# Adapted from code originally written by Blair Shevlin for "Negative affect
# influences the computations underlying food choice in bulimia nervosa".

# This version fits 3-factor model with taste delay relative to health and umami

# Packages required
required_packages <- c(
  "DEoptim", 
  "Rcpp",
  "parallel", 
  "here", 
  "fs",
  "RcppParallel",
  "stats4",
  "pracma",
  "runjags",
  "tidyverse",
  "loo",
  "coda"
)

library(readxl)
library(parallel)
# Check and install missing packages
install_if_missing <- function(p) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }
}

# Install missing packages
invisible(sapply(required_packages, install_if_missing))

# Load all packages
invisible(sapply(required_packages, library, character.only = TRUE))

# Configure RcppParallel (no outer fold-level parallelism anymore, folds run sequentially)
RcppParallel::setThreadOptions(numThreads = 1) # avoid TBB thread pool + fork() (mclapply) incompatibility

base_path = "~/Documents/Github"


scriptFolder = path(base_path) / "src"
resFolder = path(base_path) / "results"

dat <- read.csv("~/Documents/Github/mTurk_July2025.csv")
source("~/Documents/Github/src/cleanup.R")
master <- t
# recode some things

########## Choice phase ########## 
master$SubID <- as.factor(master$subjectId)

master$choice.bin <- ifelse(master$choice.rating >= 6, 1, ifelse(master$choice.rating <= 4, 0, NA))

data <- master
ref <- data[is.na(data$choice.rating), c("subjectId", "Umami_og", "Healthiness_og", "Tastiness_og")]
colnames(ref) <- c("subjectId", "ref_p", "ref_h", "ref_t")
data <- dplyr::left_join(data, ref, by = "subjectId")
data <- subset(data, !is.na(data$choice.bin))
data$choice <- ifelse(data$choice.bin==1,1,0)
data$rt <- data$choice.rt.unscale/1000
data$subjID <- data$subjectId

data$vt_r    <- data$ref_t
data$vh_r    <- data$ref_h
data$vt_l    <- data$Tastiness_og
data$vh_l    <- data$Healthiness_og
data$subjectId <- data$SubID

# Use absolute ratings (not relative to neutral)
data$tasterating <- data$vt_l - data$vt_r
data$healthrating <- data$vh_l - data$vh_r
data$umamirating <- data$Umami_og - data$ref_p

fi <- data
name <- "mTurk"
ALL <- fi %>% dplyr::select(subjectId, tasterating, healthrating, umamirating, rt, choice)
ALL[,2:5] <- sapply(ALL[,2:5], as.numeric)

ALL <- subset(ALL, !is.na(ALL$choice))

Data <- ALL

Data <- subset(Data, !is.na(Data$tasterating))
Data <- subset(Data, !is.na(Data$healthrating))
Data <- subset(Data, !is.na(Data$umamirating))

Data <- subset(Data, !(Data$rt < 0.2))
Data <- subset(Data, !(Data$rt > 5))

Data <- subset(Data, !is.na(Data$rt))
trial <- Data %>% group_by(subjectId) %>% count()
trial <- subset(trial, trial$n > 35)
Data <- subset(Data,Data$subjectId %in% trial$subjectId)

# RT is positive if left food item choosen, negative if right food item chosen
idx = which(Data$choice==0)
Data$RT <- Data$rt
Data$RT[idx] = Data$rt[idx] * -1
Data$subject <- Data$subjectId

# Scale values within subject
Data <- Data %>% group_by(subjectId) %>%  mutate(hd = scale(healthrating)[,1])
Data <- Data %>% group_by(subjectId) %>%  mutate(td = scale(tasterating)[,1])
Data <- Data %>% group_by(subjectId) %>%  mutate(pd = scale(umamirating)[,1])

idxP = as.numeric(ordered(Data$subjectId)) #makes a sequentially numbered subj index
rtpos = Data$rt

subs = Data %>% dplyr::select(subjectId) %>% unique()
write.csv(subs,"results/new_Sept_3fac_delaytaste_sublist.csv")

set.seed(42)

# Define number of folds (or train/test split)
K <- 5
Data$fold <- rep(1:K, length.out = nrow(Data))

# Simulation code (load once before parallel processing)
# Note: M6 uses a time-varying 3-factor model (taste, health, umami with taste timing)
sourceCpp("src/ddm_cpp_M6.cpp") 


# CRPS: proper scoring rule for a predictive SAMPLE vs one observation.
#   CRPS = mean|x_i - y| - (1/n^2) * sum_i (2i-n-1) * x_(i)      [x sorted]
# Equivalent to the O(n^2) double-sum form but O(n log n). Applied to SIGNED RT,
# so the sign carries the choice and the magnitude carries the RT -- unlike the
# choice-only MSE, which discards the RT the timing models differ on.
crps_sample <- function(x, y) {
  x <- sort(x[is.finite(x)]); n <- length(x)
  if (n == 0L) return(NA_real_)
  mean(abs(x - y)) - sum((2 * seq_len(n) - n - 1) * x) / (n * n)
}

# Additional scores, all computed from the SAME simulated draws (the expensive
# part is generating them, so extra metrics are effectively free).
#   crps_signed  : sign carries choice, magnitude carries RT. Implicitly prices
#                  one choice error at ~2*|RT| seconds of timing error.
#   crps_absrt   : CRPS on |RT| among draws matching the OBSERVED choice, so
#                  timing is scored without any choice/RT exchange rate.
#   brier        : squared error on the predicted choice probability (proper).
#   qdiff        : mean |predicted - observed| RT quantile gap, computed within
#                  the observed response type -- the DDM quantile-probability idea.
QPROBS <- c(.1, .3, .5, .7, .9)
score_draws <- function(sim, rt_obs) {
  ok   <- is.finite(sim)
  s    <- sim[ok]
  n    <- length(s)
  up   <- rt_obs > 0
  out  <- c(brier = NA_real_, crps_absrt = NA_real_, qdiff = NA_real_, p_match = NA_real_)
  if (n == 0L) return(out)
  p_up            <- mean(s > 0)
  out[["brier"]]  <- (as.numeric(up) - p_up)^2
  same            <- if (up) s[s > 0] else s[s < 0]      # draws on the observed side
  out[["p_match"]] <- length(same) / n
  if (length(same) >= 5L) {
    out[["crps_absrt"]] <- crps_sample(abs(same), abs(rt_obs))
    out[["qdiff"]]      <- mean(abs(quantile(abs(same), QPROBS, names = FALSE) - abs(rt_obs)))
  }
  out
}

SIM_HORIZON <- 8   # must match `double T = ...` in the .cpp; censored draws return exactly this

# Cross-validation
cat("Starting cross-validation...\n")
fold_results <- vector("list", 5)
for (k in 1:5) {
  set.seed(1000 + k)  # per-fold seed: makes each fold reproducible on its own
  rm(list=c("train_data","test_data","hd","td","pd","rtpos","ns","idxP","N","y","dat"))
  train_data <- subset(Data,! Data$fold==k) #Data %>% filter(fold != k)
  test_data  <- subset(Data, Data$fold==k)  #Data %>% filter(fold == k)
  idxP = as.numeric(ordered(train_data$subjectId))
  y=train_data$RT
  N=length(train_data$RT)
  idxP=idxP
  hd=train_data$hd
  td=train_data$td
  pd=train_data$pd
  rtpos = train_data$rt
  ns = length(unique(idxP))
  
dat <- dump.format(list(N=N, y=y, idxP=idxP, hd=hd, td=td, pd=pd, rt=rtpos, ns=ns))

inits3 <- dump.format(list( alpha.mu=2,
                            alpha.pr=0.5, theta.mu=0.1,
                            theta.pr=0.05,  b1.mu=0.3, b1.pr=0.05, b2.mu=0.05, b2.pr=0.05, b3.mu=0.01, b3.pr=0.05,
                            time_taste.mu=0.1, time_taste.pr=0.05,
                            bias.mu=0.4,
                            bias.kappa=1, y_pred=y,  .RNG.name="base::Super-Duper", .RNG.seed=99999))


inits2 <- dump.format(list( alpha.mu=2.2,
                            alpha.pr=0.05,  theta.mu=0.01,
                            theta.pr=0.05, b1.mu=0.3, b1.pr=0.05, b2.mu=0.1, b2.pr=0.05, b3.mu=0.1, b3.pr=0.05,
                            time_taste.mu=0.2, time_taste.pr=0.001,
                            bias.mu=0.4,
                            bias.kappa=1, y_pred=y,  .RNG.name="base::Wichmann-Hill", .RNG.seed=1234))

inits1 <- dump.format(list( alpha.mu=2.4,
                            alpha.pr=0.05, theta.mu=0.15,
                            theta.pr=0.05, b1.mu=0.1, b1.pr=0.05, b2.mu=0.05, b2.pr=0.05, b3.mu=0.05, b3.pr=0.05,
                            time_taste.mu=0, time_taste.pr=0.01,
                            bias.mu=0.4,
                            bias.kappa=1, y_pred=y, .RNG.name="base::Mersenne-Twister", .RNG.seed=6666 ))
      
      monitor = c(
        "alpha.mu","theta.mu",
        "b1.mu","b2.mu","b3.mu", "bias.mu",
        "time_taste.mu",
        "b1.p","b2.p","b3.p",
        "time_taste.p",
        "theta.p",
        "bias",
        "alpha.p")
      
      model = "src/ddm_priors_M6.txt"
      
     fit <- run.jags(model=model, 
                          monitor=monitor, data=dat, n.chains=3, inits=c(inits1,inits2, inits3), 
                          plots = TRUE, method="parallel", module="wiener", burnin=50000, sample=10000, thin=1)
      
      
      samples <- as.mcmc.list(fit)
      combined_samples <- do.call(rbind, lapply(samples, as.matrix))  # dims: (iterations * chains) × parameters
      idx_samp <- sample(nrow(combined_samples), 1000, replace = FALSE)
      subsampled_samples <- combined_samples[idx_samp, ]
      
      
      # Assuming parameter names are like "alpha[1]", "alpha[2]", ...
      # Get subject indices in training set (must match order used in JAGS)
      train_subjects <- unique(train_data$subject)
      
      n_iter <-1000
      
      get_param_matrix <- function(param_name, subsample = subsampled_samples) {
        param_cols <- grep(paste0("^", param_name, "\\.p\\["), colnames(subsample))
        subsample[, param_cols, drop = FALSE]  # returns matrix: 1000 draws × subjects
      }
      
      get_param_matrix_bias <- function(param_name, subsample = subsampled_samples) {
        param_cols <- param_cols <- grep("^bias\\[[0-9]+\\]$", colnames(subsample))
        subsample[, param_cols, drop = FALSE]  # returns matrix: 1000 draws × subjects
      }
      
      alpha_samples <- get_param_matrix("alpha")
      b1_samples    <- get_param_matrix("b1")
      b2_samples    <- get_param_matrix("b2")
      b3_samples    <- get_param_matrix("b3")
      theta_samples <- get_param_matrix("theta")
      bias_samples  <- get_param_matrix_bias("bias")
      time_taste_samples <- get_param_matrix("time_taste")
      
      # Map test subjects to indices in training subjects
      # If test subjects are NOT in training, we cannot get individual params!
      test_subjects <- unique(test_data$subject)
      test_in_train_idx <- match(test_subjects, train_subjects)
      
      if (any(is.na(test_in_train_idx))) {
        stop("Some test subjects are not in training data; individual-level prediction not possible.")
      }
      
      # Calculate predictive log-likelihood per test subject
      # collectors for the raw simulation draws (the expensive artefact)
      draw_store <- vector("list", 0)
      meta_store <- vector("list", 0)
      subject_diffs <- numeric(length(test_subjects))
      subject_crps  <- numeric(length(test_subjects))
      subject_pmatch <- numeric(length(test_subjects))
      subject_qdiff <- numeric(length(test_subjects))
      subject_crps_abs <- numeric(length(test_subjects))
      subject_brier <- numeric(length(test_subjects))
      subject_cens  <- numeric(length(test_subjects))
      names(subject_diffs) <- test_subjects
      names(subject_crps)  <- test_subjects
      names(subject_pmatch) <- test_subjects
      names(subject_qdiff) <- test_subjects
      names(subject_crps_abs) <- test_subjects
      names(subject_brier) <- test_subjects
      names(subject_cens)  <- test_subjects
      
      for (subj in test_subjects) {
        print(subj)
        subj_trials <- test_data %>% filter(subject == subj)
        subj_idx <- match(subj, train_subjects)
        
        trial_diffs <- numeric(nrow(subj_trials))
        trial_crps  <- numeric(nrow(subj_trials))
        trial_pmatch <- numeric(nrow(subj_trials))
        trial_qdiff <- numeric(nrow(subj_trials))
        trial_crps_abs <- numeric(nrow(subj_trials))
        trial_brier <- numeric(nrow(subj_trials))
        trial_cens  <- numeric(nrow(subj_trials))
        
        for (i in seq_len(nrow(subj_trials))) {
          rt_i <- subj_trials$RT[i]
          
          subframe <- data.frame(cbind(alpha_samples[,subj_idx],b1_samples[,subj_idx],b2_samples[,subj_idx],b3_samples[,subj_idx],theta_samples[,subj_idx],bias_samples[,subj_idx],time_taste_samples[,subj_idx]))
          colnames(subframe) <- c("alpha","b1","b2","b3","theta","bias","time_taste")
          hd <- subj_trials$hd[i]
          td <- subj_trials$td[i]
          pd <- subj_trials$pd[i]
          
          # Extract all posterior samples of parameters for this subject
          sim_samples <- numeric(n_iter)
          seeds <- sample.int(2147483647L, n_iter)  # reproducible seeds for the C++ simulator
          
          for (iter in 1:n_iter) {
            b1=subframe$b1[iter]
            b2=subframe$b2[iter]
            b3=subframe$b3[iter]
            alpha=subframe$alpha[iter]
            theta=subframe$theta[iter]
            bias=subframe$bias[iter]
            tIn_t=subframe$time_taste[iter]
            sim_samples[iter] <- ddm3t_parallel(d_v = b1,d_h = b2,d_p = b3,thres = alpha,nDT = theta,bias = bias,vd =td ,hd =hd,pd =pd,tIn_t=tIn_t,N=1,sd_n=1,seed = seeds[iter])
          }
          
          y_no <- ifelse(sim_samples > 0,1,0)
          y_no_true <- ifelse(rt_i > 0,1,0)
          mse <- mean((y_no_true - y_no)^2)
          trial_diffs[i] <- mse
          trial_crps[i]  <- crps_sample(sim_samples, rt_i)
          trial_cens[i]  <- mean(sim_samples == SIM_HORIZON)
          .sc <- score_draws(sim_samples, rt_i)
          trial_brier[i]    <- .sc[["brier"]]
          trial_crps_abs[i] <- .sc[["crps_absrt"]]
          trial_qdiff[i]    <- .sc[["qdiff"]]
          trial_pmatch[i]   <- .sc[["p_match"]]
          draw_store[[length(draw_store) + 1L]] <- sim_samples
          meta_store[[length(meta_store) + 1L]] <- data.frame(
            subjectId = as.character(subj), trial = i,
            td = td, hd = hd, RT_obs = rt_i, choice_obs = as.integer(rt_i > 0))
        }
        
        subject_diffs[subj] <- mean(trial_diffs)
        subject_crps[subj]  <- mean(trial_crps, na.rm = TRUE)
        subject_cens[subj]  <- mean(trial_cens)
        subject_brier[subj]    <- mean(trial_brier,    na.rm = TRUE)
        subject_crps_abs[subj] <- mean(trial_crps_abs, na.rm = TRUE)
        subject_qdiff[subj]    <- mean(trial_qdiff,    na.rm = TRUE)
        subject_pmatch[subj]   <- mean(trial_pmatch,   na.rm = TRUE)
      }
      
      # Return results for this fold
      fold_result <- list(
        test_subject = subject_diffs,
        crps         = subject_crps,
        censored     = subject_cens,
        brier        = subject_brier,
        crps_abs     = subject_crps_abs,
        qdiff        = subject_qdiff,
        pmatch       = subject_pmatch
      )
      
      # Save individual fold result
      a <- data.frame(fold_result)
      a$subid <- rownames(a)
      write.csv(a, paste0("results/new_Sept_3fac_delaytaste_cross_validated_fold_", k, ".csv"))

      # --- save raw draws + metrics for this fold (lets any future metric be
      # --- recomputed without re-running the simulation) ---
      sim_mat <- do.call(rbind, draw_store)          # n_trials x n_iter
      saveRDS(list(model = "M6", fold = k,
                   sim_draws = sim_mat,
                   trials    = do.call(rbind, meta_store),
                   metrics   = fold_result,
                   horizon   = SIM_HORIZON, n_iter = n_iter),
              file = sprintf("results/fits/cv_M6_fold%d.rds", k))
      cat(sprintf("  fold %d: saved %d x %d draw matrix (%.0f MB)\n", k,
                  nrow(sim_mat), ncol(sim_mat), object.size(sim_mat)/1048576))
      rm(sim_mat, draw_store, meta_store); invisible(gc())
      fold_results[[k]] <- fold_result
} # End of fold loop

cat("Parallel cross-validation completed!\n")


s <- subs
for(i in 1:5){
  a <- data.frame(fold_results[[i]])
  a$subid <- rownames(a)
names(a)[names(a) != "subid"] <- paste0(names(a)[names(a) != "subid"], "_f", i)  # keep folds distinguishable
  s <- merge(s,a,by.x=1,by.y="subid")
}

write.csv(s,"results/new_Sept_3fac_delaytaste_cross_validated_folds.csv")
