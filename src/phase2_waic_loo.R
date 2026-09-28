# Phase 2: marginal WAIC / PSIS-LOO for M1-M6.
#
# The log_lik JAGS reports is conditional on the per-trial latent drift w[i],
# which is fitted to that same trial -- so LOO/WAIC on it are degenerate.
# Here we integrate w[i] out by Gauss-Hermite quadrature (see
# new/wiener_marginal.cpp, validated against RWiener to ~1e-15) to get a
# genuine pointwise predictive density, then feed that to loo.
#
# li.hat is reconstructed per model to match src/ddm_priors_M*.txt exactly.
# NOTE: JAGS step(x) == 1 when x >= 0, else 0.

suppressMessages({library(Rcpp); library(statmod); library(loo)})
setwd("/home/processor/Documents/Github")
sourceCpp("new/wiener_marginal.cpp")

NQ     <- 25                      # quadrature nodes
NDRAWS <- 3000                    # posterior draws to use (all of them)
gq     <- statmod::gauss.quad(NQ, "hermite")

# which time parameter each model carries (NA = no time-varying term)
time_par <- c(M1 = NA, M2 = NA, M3 = "time.p",
              M4 = "time_umami.p", M5 = "time_health.p", M6 = "time_taste.p")
has_b3   <- c(M1 = FALSE, M2 = TRUE, M3 = FALSE, M4 = TRUE, M5 = TRUE, M6 = TRUE)

loglik_matrix <- function(m) {
  fit  <- readRDS(sprintf("results/fits/fullfit_M%d.rds", m))
  D    <- fit$Data
  post <- fit$post
  mn   <- paste0("M", m)

  # subject index exactly as JAGS built it
  idx <- as.numeric(ordered(D$subjectId))
  ci  <- function(nm) {
    j <- match(paste0(nm, "[", idx, "]"), colnames(post))
    if (any(is.na(j))) stop(sprintf("%s: cannot locate columns for %s", mn, nm))
    j
  }
  c_alpha <- ci("alpha.p"); c_theta <- ci("theta.p"); c_bias <- ci("bias")
  c_etau  <- ci("e.p.tau"); c_b1 <- ci("b1.p"); c_b2 <- ci("b2.p")
  c_b3    <- if (has_b3[[mn]]) ci("b3.p") else NULL
  c_time  <- if (!is.na(time_par[[mn]])) ci(time_par[[mn]]) else NULL

  S  <- min(NDRAWS, nrow(post))
  sel <- round(seq(1, nrow(post), length.out = S))
  n  <- nrow(D)
  ll <- matrix(NA_real_, nrow = S, ncol = n)

  tpos  <- D$rt                    # positive RT (JAGS 'rt', used for dT)
  upper <- D$RT > 0                # signed RT carries the choice
  td <- D$td; hd <- D$hd
  pd <- if (has_b3[[mn]]) D$pd else NULL

  for (i in seq_len(S)) {
    s  <- sel[i]
    a  <- post[s, c_alpha]; ta <- post[s, c_theta]; be <- post[s, c_bias]
    sdv <- 1 / sqrt(post[s, c_etau])

    vt <- post[s, c_b1] * td
    vh <- post[s, c_b2] * hd
    vp <- if (has_b3[[mn]]) post[s, c_b3] * pd else rep(0, n)

    if (is.na(time_par[[mn]])) {
      mu <- vt + vh + vp                                   # M1, M2: static
    } else {
      tm      <- post[s, c_time]
      dT      <- tpos - ta
      tstep   <- (dT - abs(tm)) / dT
      smaller <- as.numeric(tstep >= 0)                    # JAGS step()
      f       <- as.numeric(tm >= 0)                        # JAGS step()
      disc    <- smaller * tstep                            # discount applied to the delayed attribute

      if (mn == "M3") {          # f=1 -> health delayed;  f=0 -> taste delayed
        mu <- (f * vt + (1 - f) * disc * vt) +
              ((1 - f) * vh + f * disc * vh)
      } else if (mn == "M4") {   # f=1 -> taste+health delayed; f=0 -> umami delayed
        mu <- ((1 - f) * vt + f * disc * vt) +
              ((1 - f) * vh + f * disc * vh) +
              (f * vp + (1 - f) * disc * vp)
      } else if (mn == "M5") {   # f=1 -> health delayed
        mu <- (f * vt + (1 - f) * disc * vt) +
              ((1 - f) * vh + f * disc * vh) +
              (f * vp + (1 - f) * disc * vp)
      } else if (mn == "M6") {   # f=1 -> taste delayed
        mu <- ((1 - f) * vt + f * disc * vt) +
              (f * vh + (1 - f) * disc * vh) +
              (f * vp + (1 - f) * disc * vp)
      }
    }

    ll[i, ] <- wiener_marginal_loglik(tpos, upper, a, ta, be, mu, sdv,
                                      gq$nodes, gq$weights)
  }
  list(ll = ll, n = n, S = S, model = mn)
}

res <- list()
for (m in 1:6) {
  t0 <- Sys.time()
  o  <- loglik_matrix(m)
  ll <- o$ll

  bad <- !is.finite(ll)
  cat(sprintf("\n%s: %d draws x %d trials | non-finite loglik cells: %d (%.4f%%)\n",
              o$model, o$S, o$n, sum(bad), 100 * mean(bad)))

  # trials where the density underflows for some draws are unusable for LOO
  drop <- which(colSums(bad) > 0)
  if (length(drop)) {
    cat(sprintf("  dropping %d trial(s) with any non-finite value\n", length(drop)))
    ll <- ll[, -drop, drop = FALSE]
  }

  # chains were rbind-ed in contiguous blocks, so recover chain id for r_eff
  cid  <- rep(1:3, each = ceiling(o$S / 3))[seq_len(o$S)]
  reff <- relative_eff(exp(ll), chain_id = cid)

  w <- suppressWarnings(waic(ll))
  l <- suppressWarnings(loo(ll, r_eff = reff))

  res[[o$model]] <- list(waic = w, loo = l, n_used = ncol(ll), dropped = length(drop))
  cat(sprintf("  elapsed %.1f s | elpd_waic %.1f | elpd_loo %.1f | p_loo %.1f | bad k>0.7: %d\n",
              as.numeric(difftime(Sys.time(), t0, units = "secs")),
              w$estimates["elpd_waic","Estimate"], l$estimates["elpd_loo","Estimate"],
              l$estimates["p_loo","Estimate"], sum(l$diagnostics$pareto_k > 0.7)))
  rm(o, ll); invisible(gc())
}

saveRDS(res, "results/fits/phase2_waic_loo.rds")

cat("\n\n================ MODEL COMPARISON (marginal LOO) ================\n")
tab <- data.frame(
  model    = names(res),
  n_trials = sapply(res, function(x) x$n_used),
  elpd_loo = sapply(res, function(x) x$loo$estimates["elpd_loo","Estimate"]),
  se       = sapply(res, function(x) x$loo$estimates["elpd_loo","SE"]),
  p_loo    = sapply(res, function(x) x$loo$estimates["p_loo","Estimate"]),
  elpd_waic= sapply(res, function(x) x$waic$estimates["elpd_waic","Estimate"]),
  bad_k    = sapply(res, function(x) sum(x$loo$diagnostics$pareto_k > 0.7)),
  row.names = NULL)
tab <- tab[order(-tab$elpd_loo), ]
print(tab, row.names = FALSE, digits = 6)

cat("\n--- pairwise vs best (loo_compare; only valid where trial sets match) ---\n")
nsame <- sapply(res, function(x) x$n_used)
if (length(unique(nsame)) == 1) {
  print(loo_compare(lapply(res, function(x) x$loo)))
} else {
  cat("  trial sets differ across models (different numbers dropped):\n")
  print(nsame)
  cat("  -> loo_compare not directly valid; recompute on the common trial set if needed.\n")
}
