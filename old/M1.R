# fit_stddm_models.R - Code for fitting models from "Negative affect influences the computations underlying food choice in bulimia nervosa" 
#
# Copyright (C) 2024 Blair Shevlin, <blairshevlin@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# Record of Revisions
#
# Date            Programmers                         Descriptions of Change
# ====         ================                       ======================
# 06/18/24      Blair Shevlin                         Wrote code for original manuscript


# This version fits foodType WITHIN a condition

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
  "coda",
  "foreach",
  "doParallel"
)

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

# Set up parallel processing
n_cores <- detectCores() - 1  # Leave one core free
cat("Using", n_cores, "cores for parallel processing\n")

# Set up parallel backend for foreach
registerDoParallel(cores = n_cores)

# Configure RcppParallel to use fewer threads per process since we're parallelizing folds
RcppParallel::setThreadOptions(numThreads = max(1, floor(detectCores() / 5))) # Divide cores among folds


# === DATA LOADING ===
# Replace these paths with your own data files
# dat <- read.csv("path/to/your/data.csv")
# source("path/to/your/cleanup_script.R")
# master <- t
# recode some things

########## Choice phase ########## 
master$SubID <- as.factor(master$subjectId)

#master$fat.e <- as.factor(master$fat.e)

master$choice.bin <- ifelse(master$choice.rating >= 6, 1, ifelse(master$choice.rating <= 4, 0, NA))
#master$choice.bin <- ifelse(master$choice.rating > 5 & master$choice.rating < 6,NA,master$choice.bin)


data <- master
data$ref_p <- data[is.na(data$choice.rating),]$Umami_og
data$ref_h <- data[is.na(data$choice.rating),]$Healthiness_og
data$ref_t <- data[is.na(data$choice.rating),]$Tastiness_og
data <- subset(data, !is.na(data$choice.bin))
data$choice <- ifelse(data$choice.bin==1,1,0)
data$rt <- data$choice.rt.unscale/1000
data$subjID <- data$subjectId

data$vt_r    <- data$ref_t
data$vh_r    <- data$ref_h
data$vt_l    <- data$Tastiness_og
data$vh_l    <- data$Healthiness_og
data$subjectId <- data$SubID


#data$tasterating <- data$vt_l - data$vt_r
#data$healthrating <- data$vh_l - data$vh_r

#do not need to standardize within individual OR compare with neutral item for difference since
# the neutral item is constant
data$tasterating <- data$vt_l - data$vt_r
data$healthrating <- data$vh_l - data$vh_r
data$umamirating <- data$Umami_og - data$ref_p

fi <- data
name <- "mTurk"
ALL <- fi %>% dplyr::select(subjectId, tasterating, healthrating, rt, choice)
ALL[,2:5] <- sapply(ALL[,2:5], as.numeric)

ALL <- subset(ALL, !is.na(ALL$choice))
ALL <- subset(ALL, !is.na(ALL$rt))


Data <- ALL

#Data$td     <- Data$vt_l - Data$vt_r
#Data$hd     <- Data$vh_l - Data$vh_r
Data <- subset(Data, !is.na(Data$tasterating))
Data <- subset(Data, !is.na(Data$healthrating))

Data <- subset(Data, !(Data$rt < 0.2))
Data <- subset(Data, !(Data$rt > 5))

Data <- subset(Data, !is.na(Data$rt))
trial <- Data %>% group_by(subjectId) %>% count()
trial <- subset(trial, trial$n > 35)
Data <- subset(Data,Data$subjectId %in% trial$subjectId)

#Data = Data[Data$rt>0.2,]

# RT is positive if left food item choosen, negative if right food item chosen
idx = which(Data$choice==0)

Data$RT <- Data$rt
Data$RT[idx] = Data$rt[idx] * -1
Data$subject <- Data$subjectId

# scale value of health and taste difference
#Data$hds = scale(Data$hd)[,1]
#Data$tds = scale(Data$td)[,1]

#hd = scale(Data$healthrating,center=F)[,1]
#td = scale(Data$tasterating,center=F)[,1]

Data <- Data %>% group_by(subjectId) %>%  mutate(hd = scale(healthrating)[,1])
Data <- Data %>% group_by(subjectId) %>%  mutate(td = scale(tasterating)[,1])

hd = Data$healthrating
td = Data$tasterating

idxP = as.numeric(ordered(Data$subjectId)) #makes a sequentially numbered subj index
rtpos = Data$rt

subs = Data %>% dplyr::select(subjectId) %>% unique()
write.csv(subs,"results/Sept_noTV_sublist.csv")

y= Data$RT
N = length(y)
ns = length(unique(idxP))        


set.seed(42)


# Define number of folds (or train/test split)
K <- 5
Data$fold <- rep(1:K, length.out = nrow(Data))






# Fit initial
fold_results <- NULL
# Simulation code (load once before parallel processing)
sourceCpp("src/ddm_cpp_M1.cpp")

# Parallel cross-validation
cat("Starting parallel cross-validation with", n_cores, "cores...\n")
fold_results <- foreach(k = 1:5, .combine = list, .multicombine = TRUE, .packages = c("runjags", "R2jags", "coda")) %dopar% {
  train_data <- subset(Data,! Data$fold==k) #Data %>% filter(fold != k)
  test_data  <- subset(Data, Data$fold==k)  #Data %>% filter(fold == k)
  idxP = as.numeric(ordered(train_data$subjectId)) 
  y=train_data$RT
  N=length(train_data$RT)
  idxP=idxP
  hd=train_data$hd
  td=train_data$td
  rtpos = train_data$rt
  ns = length(unique(idxP))
  
  
  jags_data <- dump.format(list(N=N, y=y, idxP=idxP, hd=hd, td=td, rt=rtpos, ns=ns))
      
      inits3 <- dump.format(list( alpha.mu=2,  
                                  alpha.pr=0.5, theta.mu=0.1,
                                  theta.pr=0.05,  b1.mu=0.3, b1.pr=0.05, b2.mu=0.01, b2.pr=0.05, 
                                  bias.mu=0.4,
                                  bias.kappa=1, y_pred=y,  .RNG.name="base::Super-Duper", .RNG.seed=99999))
      
      
      inits2 <- dump.format(list( alpha.mu=2.2, 
                                  alpha.pr=0.05, theta.mu=0.01,
                                  theta.pr=0.05, b1.mu=0.3, b1.pr=0.05, b2.mu=0.1, b2.pr=0.05, 
                                  bias.mu=0.4,
                                  bias.kappa=1, y_pred=y,  .RNG.name="base::Wichmann-Hill", .RNG.seed=1234))
      
      inits1 <- dump.format(list( alpha.mu=2.4,
                                  alpha.pr=0.05, theta.mu=0.15,
                                  theta.pr=0.05, b1.mu=0.1, b1.pr=0.05, b2.mu=0.05, b2.pr=0.05, 
                                  bias.mu=0.4,
                                  bias.kappa=1, y_pred=y, .RNG.name="base::Mersenne-Twister", .RNG.seed=6666 ))
      
      monitor = c(
        "alpha.mu","theta.mu",
        "b1.mu","b2.mu","bias.mu",
        "b1.p","b2.p", 
        "theta.p", 
        "bias",
        "alpha.p","log_lik",
        "deviance")
      
      model = "src/ddm_priors_M1.txt"
      
      fit <- runjags::run.jags(model=model, 
                                   monitor=monitor, data=jags_data, n.chains=3, inits=c(inits1,inits2, inits3), 
                                   plots = TRUE, method="parallel", module="wiener", burnin=50000, sample=10000, thin=1)
      
      
      samples <- as.mcmc.list(fit)
      samples_array <- coda::as.array.mcmc.list(samples)
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
      theta_samples <- get_param_matrix("theta")
      bias_samples  <- get_param_matrix_bias("bias")  # double-check name if not "bias.p[...]"
      
      # Map test subjects to indices in training subjects
      # If test subjects are NOT in training, we cannot get individual params!
      test_subjects <- unique(test_data$subject)
      test_in_train_idx <- match(test_subjects, train_subjects)
      
      if (any(is.na(test_in_train_idx))) {
        stop("Some test subjects are not in training data; individual-level prediction not possible.")
      }
      
      # Calculate predictive log-likelihood per test subject
      subject_diffs <- numeric(length(test_subjects))
      names(subject_diffs) <- test_subjects
      
      for (subj in test_subjects) {
        print(subj)
        subj_trials <- test_data %>% filter(subject == subj)
        subj_idx <- match(subj, train_subjects)
        
        trial_diffs <- numeric(nrow(subj_trials))
        
        for (i in seq_len(nrow(subj_trials))) {
          rt_i <- subj_trials$RT[i]
          
          subframe <- data.frame(cbind(alpha_samples[,i],b1_samples[,i],b2_samples[,i],theta_samples[,i],bias_samples[,i]))
          colnames(subframe) <- c("alpha","b1","b2","theta","bias")
          hd <- subj_trials$hd[i]
          td <- subj_trials$td[i]
          pd <- subj_trials$pd[i]
          
          # Extract all posterior samples of parameters for this subject
          sim_samples <- numeric(n_iter)
          
          for (iter in 1:n_iter) {
            b1=subframe$b1[i]
            b2=subframe$b2[i]
            alpha=subframe$alpha[i]
            theta=subframe$theta[i]
            bias=subframe$bias[i]
            sim_samples[iter] <- ddm2_parallel(d_v = b1,d_h = b2,thres = alpha,nDT = theta,bias = bias,vd =td ,hd =hd ,N=1,sd_n=1)
          }
          
          y_no <- ifelse(sim_samples > 0,1,0)
          y_no_true <- ifelse(rt_i > 0,1,0)
          mse <- mean((y_no_true - y_no)^2)
          trial_diffs[i] <- mse
        }
        
        subject_diffs[subj] <- mean(trial_diffs)
      }
 
      
      # Return results for this fold
      fold_result <- list(
        #model_fit = fit
        test_subject = subject_diffs
      )
      
      return(fold_result)
} # End of parallel foreach

cat("Parallel cross-validation completed!\n")

# Stop parallel cluster
stopImplicitCluster()

s <- subs
for(i in 1:5){
a <- data.frame(fold_results[[i]])
a$subid <- rownames(a)
s <- merge(s,a,by.x=1,by.y=2)
}

write.csv(s,"results/Sept_noTV_cross_validated_folds.csv")
