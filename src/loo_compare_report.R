# Model comparison from the marginal LOO objects, with and without the trials
# flagged by high Pareto k.
#
# Requires results/fits/phase2_waic_loo.rds, produced by src/phase2_waic_loo.R
# (which builds the marginal log_lik via src/wiener_marginal.cpp and calls loo).
#
# Why the second table: PSIS-LOO's importance sampling degrades for observations
# with pareto_k > 0.7, so those pointwise elpd values are unreliable. Rather than
# drop them per model (which would compare different trial sets), we restrict to
# the trials reliable in EVERY model, so all six are scored on identical data.

suppressMessages(library(loo))
setwd("/home/processor/Documents/Github")

res  <- readRDS("results/fits/phase2_waic_loo.rds")
mods <- names(res)

K  <- sapply(res, function(x) x$loo$diagnostics$pareto_k)      # trials x models
PW <- sapply(res, function(x) x$loo$pointwise[, "elpd_loo"])   # trials x models
stopifnot(identical(dim(K), dim(PW)))

good <- apply(K, 1, function(k) all(k <= 0.7))
cat(sprintf("trials: %d | reliable (k <= 0.7) in ALL models: %d (%.1f%%)\n\n",
            nrow(K), sum(good), 100 * mean(good)))

full <- colSums(PW)
sub  <- colSums(PW[good, , drop = FALSE])
se_f <- apply(PW, 2, function(v) sqrt(length(v)) * sd(v))
se_s <- apply(PW[good, , drop = FALSE], 2, function(v) sqrt(length(v)) * sd(v))

tab <- data.frame(
  model      = mods,
  elpd_loo   = full,
  se         = se_f,
  d_best     = full - max(full),
  bad_k      = colSums(K > 0.7),
  elpd_good  = sub,
  se_good    = se_s,
  d_best_good = sub - max(sub),
  row.names  = NULL)
tab <- tab[order(-tab$elpd_loo), ]

cat("========== ALL TRIALS ==========\n")
print(tab[, c("model","elpd_loo","se","d_best","bad_k")], row.names = FALSE, digits = 6)
cat("\n========== HIGH-PARETO-K TRIALS EXCLUDED ==========\n")
print(tab[, c("model","elpd_good","se_good","d_best_good")], row.names = FALSE, digits = 6)

cat("\nranking (all trials)  :", paste(mods[order(-full)], collapse = " > "), "\n")
cat("ranking (k <= 0.7)    :", paste(mods[order(-sub)],  collapse = " > "), "\n")
cat("rankings agree        :", identical(order(-full), order(-sub)), "\n")

# paired SE on the difference vs the best model -- the standard loo_compare
# quantity, but recomputed here so it can also be done on the restricted set
best_f <- mods[which.max(full)]; best_s <- mods[which.max(sub)]
cat(sprintf("\npairwise vs %s (paired SE of the pointwise difference):\n", best_f))
cat(sprintf("  %-5s %12s %9s %8s | %12s %9s\n",
            "model","d_elpd_all","se_diff","z","d_elpd_good","se_diff"))
for (m in mods[order(-full)]) {
  if (m == best_f) next
  d_all  <- PW[, m] - PW[, best_f]
  d_good <- PW[good, m] - PW[good, best_s]
  se_all  <- sqrt(length(d_all))  * sd(d_all)
  se_good <- sqrt(length(d_good)) * sd(d_good)
  cat(sprintf("  %-5s %12.1f %9.2f %8.1f | %12.1f %9.2f\n",
              m, sum(d_all), se_all, sum(d_all)/se_all, sum(d_good), se_good))
}

cat("\nloo_compare() on the full objects, for reference:\n")
print(loo_compare(lapply(res, function(x) x$loo)))

write.csv(tab, "results/fits/loo_comparison_table.csv", row.names = FALSE)
cat("\nwritten: results/fits/loo_comparison_table.csv\n")
