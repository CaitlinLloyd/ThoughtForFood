# Parameter recovery for M2 (the winning model: 3-factor, no time-varying term).
#
# Design decisions, and why:
#  * Generating values are the per-subject POSTERIOR MEANS from the full-data fit,
#    and recovery is judged against the refit's posterior means -- like for like.
#    This matches the BN reference implementation.
#    Caveat to carry into interpretation: posterior means are shrunk toward the
#    group mean (measured on this fit, they retain 53% of the between-subject
#    spread for b3, 69% for b2, but ~99% for theta and 97% for bias). So recovery
#    is being tested over a compressed range for the weight parameters, and will
#    look better there than it would for the full population spread.
#  * The simulator applies that drift directly with no trial-level drift noise --
#    matching the reference implementation, which calls ddm2_parallel with point
#    values and sd_n=1. So we recover under a slightly simplified generator, and
#    e.p.tau is not part of the recovery test.
#  * Simulated RTs are rejected outside 0.2 < |rt| < 5 to match the filter applied
#    to the REAL data. So this tests the pipeline as actually run, filter included --
#    the JAGS model has no truncation term, so any systematic offset is a property
#    of the pipeline and would be present in the real estimates too.
#  * The rejection loop uses repeat/break (upstream's `while` retains an unchecked
#    redraw) and has an attempt cap so a pathological parameter set cannot spin.
#  * numThreads=1 is load-bearing: the worker holds one shared RNG, so N>1 across
#    threads is a data race.

suppressMessages({library(Rcpp); library(RcppParallel); library(runjags); library(dplyr)
                  library(coda)})   # as.mcmc.list lives here -- omitting it killed a completed 3h run
setThreadOptions(numThreads = 1)
setwd("/home/processor/Documents/Github")
sourceCpp("src/ddm_cpp_M2.cpp")

set.seed(20250819)
RT_LO <- 0.2; RT_HI <- 5; MAXATT <- 200
HORIZON <- 8            # matches `double T = 8` in the .cpp

# ---------------------------------------------------------------- generating values
fit0 <- readRDS("results/fits/fullfit_M2.rds")
D    <- fit0$Data
post <- fit0$post
idx  <- as.numeric(ordered(D$subjectId)); ns <- max(idx)

# Generating values = per-subject posterior MEANS, compared against the refit's
# posterior means (like for like, matching the BN reference). Note these are
# shrunk toward the group mean, so the weight parameters are recovered over a
# compressed range -- see the header note.
pm <- function(nm) colMeans(post[, paste0(nm, "[", 1:ns, "]")])

truth <- data.frame(subj = 1:ns,
                    alpha = pm("alpha.p"), theta = pm("theta.p"), bias = pm("bias"),
                    b1 = pm("b1.p"), b2 = pm("b2.p"), b3 = pm("b3.p"))
cat(sprintf("generating values = per-subject posterior means from %d draws\n", nrow(post)))
cat(sprintf("generating values from %d subjects, %d trials\n", ns, nrow(D)))
print(round(sapply(truth[,-1], function(x) c(min=min(x), med=median(x), max=max(x))), 3))

# ---------------------------------------------------------------- simulate
sim_rt <- numeric(nrow(D)); n_cens <- 0L; n_reject <- 0L; n_fail <- 0L
for (i in seq_len(nrow(D))) {
  p <- idx[i]
  for (att in seq_len(MAXATT)) {
    v <- ddm2_parallel(d_v = truth$b1[p], d_h = truth$b2[p], d_p = truth$b3[p],
                       thres = truth$alpha[p], nDT = truth$theta[p], bias = truth$bias[p],
                       vd = D$td[i], hd = D$hd[i], pd = D$pd[i],
                       sd_n = 1, N = 1, seed = sample.int(2147483647L, 1))
    if (v == HORIZON) { n_cens <- n_cens + 1L; n_reject <- n_reject + 1L; next }
    if (abs(v) > RT_LO && abs(v) < RT_HI) break
    n_reject <- n_reject + 1L
    if (att == MAXATT) { v <- NA_real_; n_fail <- n_fail + 1L }
  }
  sim_rt[i] <- v
}
cat(sprintf("\nsimulated %d trials | rejected %d draws (%.1f%% of attempts) | censored %d | failed %d\n",
            nrow(D), n_reject, 100*n_reject/(nrow(D)+n_reject), n_cens, n_fail))

sim <- D
sim$RT     <- sim_rt
sim$rt     <- abs(sim_rt)
sim$choice <- as.integer(sim_rt > 0)
sim <- sim[is.finite(sim$RT), ]
cat(sprintf("kept %d trials | mean|RT| %.3f (real %.3f) | P(yes) %.3f (real %.3f)\n",
            nrow(sim), mean(sim$rt), mean(D$rt), mean(sim$choice), mean(D$choice)))
saveRDS(list(truth = truth, sim = sim, n_cens = n_cens, n_reject = n_reject),
        "results/fits/recover_M2_simdata.rds")
write.csv(sim, "results/fits/recover_M2_simdata.csv", row.names = FALSE)

# ---------------------------------------------------------------- refit
idxP <- as.numeric(ordered(sim$subjectId))
y <- sim$RT; N <- length(y); ns2 <- length(unique(idxP))
dat <- dump.format(list(N=N, y=y, idxP=idxP, hd=sim$hd, td=sim$td, pd=sim$pd,
                        rt=sim$rt, ns=ns2))
inits3 <- dump.format(list(alpha.mu=2, alpha.pr=0.5, theta.mu=0.1, theta.pr=0.05,
  b1.mu=0.3,b1.pr=0.05,b2.mu=0.01,b2.pr=0.05,b3.mu=0.2,b3.pr=0.05, bias.mu=0.4,
  bias.kappa=1, y_pred=y, .RNG.name="base::Super-Duper", .RNG.seed=99999))
inits2 <- dump.format(list(alpha.mu=2.2, alpha.pr=0.05, theta.mu=0.01, theta.pr=0.05,
  b1.mu=0.3,b1.pr=0.05,b2.mu=0.1,b2.pr=0.05,b3.mu=0.2,b3.pr=0.05, bias.mu=0.4,
  bias.kappa=1, y_pred=y, .RNG.name="base::Wichmann-Hill", .RNG.seed=1234))
inits1 <- dump.format(list(alpha.mu=2.4, alpha.pr=0.05, theta.mu=0.15, theta.pr=0.05,
  b1.mu=0.1,b1.pr=0.05,b2.mu=0.05,b2.pr=0.05,b3.mu=0.2,b3.pr=0.05, bias.mu=0.4,
  bias.kappa=1, y_pred=y, .RNG.name="base::Mersenne-Twister", .RNG.seed=6666))
monitor <- c("alpha.mu","theta.mu","b1.mu","b2.mu","b3.mu","bias.mu",
             "b1.p","b2.p","b3.p","theta.p","bias","alpha.p","e.p.tau")

cat("\nrefitting on simulated data (same settings as the real fit)...\n")
fit <- run.jags(model="src/ddm_priors_M2.txt", monitor=monitor, data=dat, n.chains=3,
                inits=c(inits1,inits2,inits3), plots=FALSE, method="parallel",
                module="wiener", burnin=50000, sample=10000, thin=1)

# Extract and persist the posterior IMMEDIATELY, before any analysis code can
# fail. A previous run completed the full MCMC and then lost it on a missing
# function at exactly this point; the raw fit is dumped as a fallback so the
# expensive part is never at the mercy of downstream code.
post2 <- tryCatch({
  samples <- coda::as.mcmc.list(fit)
  comb <- do.call(rbind, lapply(samples, as.matrix))
  comb[round(seq(1, nrow(comb), length.out = min(3000, nrow(comb)))), , drop = FALSE]
}, error = function(e) {
  saveRDS(fit, "results/fits/recover_M2_RAWFIT_fallback.rds")
  stop("posterior extraction failed (raw fit dumped to recover_M2_RAWFIT_fallback.rds): ",
       conditionMessage(e))
})
saveRDS(list(post = post2, truth = truth, sim = sim), "results/fits/recover_M2_fit.rds")
cat(sprintf("posterior saved: %d draws x %d params\n", nrow(post2), ncol(post2)))

# ---------------------------------------------------------------- compare
pm2 <- function(nm) colMeans(post2[, paste0(nm, "[", 1:ns2, "]")])
rec <- data.frame(subj = 1:ns2, alpha = pm2("alpha.p"), theta = pm2("theta.p"),
                  bias = pm2("bias"), b1 = pm2("b1.p"), b2 = pm2("b2.p"), b3 = pm2("b3.p"))

# subjects kept after filtering, in the same ordering the refit used
keep_subj <- levels(ordered(sim$subjectId))
map <- match(keep_subj, levels(ordered(D$subjectId)))
tru <- truth[map, ]

out <- data.frame()
for (p in c("alpha","theta","bias","b1","b2","b3")) {
  t_ <- tru[[p]]; r_ <- rec[[p]]
  out <- rbind(out, data.frame(param = p, r = cor(t_, r_),
                               bias = mean(r_ - t_), rmse = sqrt(mean((r_ - t_)^2)),
                               true_sd = sd(t_), rec_sd = sd(r_)))
}
cat("\n================ PARAMETER RECOVERY (M2) ================\n")
print(out, row.names = FALSE, digits = 4)
saveRDS(out, "results/fits/recover_M2_summary.rds")
write.csv(cbind(subj = tru$subj, setNames(tru[,-1], paste0("true_", names(tru)[-1])),
                setNames(rec[,-1], paste0("rec_",  names(rec)[-1]))),
          "results/fits/recover_M2_params.csv", row.names = FALSE)
cat("\nsaved: recover_M2_simdata.rds/.csv, recover_M2_fit.rds, recover_M2_params.csv\n")
