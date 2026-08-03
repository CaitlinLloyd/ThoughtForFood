
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
data <- read.csv("path/to/your/data.csv")


########## Code the choice phase variables ########## 




data$choice.bin <- ifelse(data$choice.rating >= 6, 1, ifelse(data$choice.rating <= 4, 0, NA))
#data$choice.bin <- ifelse(data$choice.rating > 5 & data$choice.rating < 6,NA,data$choice.bin)



data$ref_u <- data[is.na(data$choice.rating),]$Umami_og
data$ref_h <- data[is.na(data$choice.rating),]$Healthiness_og
data$ref_t <- data[is.na(data$choice.rating),]$Tastiness_og
data <- subset(data, !is.na(data$choice.bin))

data$choice <- ifelse(data$choice.bin==1,1,0)
data$rt <- data$choice.rt.unscale/1000


data$vt_r    <- data$ref_t
data$vh_r    <- data$ref_h
data$vu_r    <- data$ref_u
data$vt_l    <- data$Tastiness_og
data$vh_l    <- data$Healthiness_og
data$vu_l   <- data$Umami_og

data$tasterating <- data$vt_l - data$vt_r
data$healthrating <- data$vh_l - data$vh_r
data$umamirating <- data$vu_l - data$vu_r

fi <- data
name <- "mTurk"
ALL <- fi %>% dplyr::select(subjectId, tasterating, healthrating,umamirating, rt, choice)
ALL[,2:5] <- sapply(ALL[,2:5], as.numeric)

ALL <- subset(ALL, !is.na(ALL$choice))
ALL <- subset(ALL, !is.na(ALL$rt))


Data <- ALL

Data <- subset(Data, !is.na(Data$tasterating))
Data <- subset(Data, !is.na(Data$healthrating))
Data <- subset(Data, !is.na(Data$tasterating))
Data <- subset(Data, !is.na(Data$umamirating))


Data <- subset(Data, !(Data$rt < 0.2))
Data <- subset(Data, !(Data$rt > 5))

Data <- subset(Data, !is.na(Data$rt))
trial <- Data %>% group_by(subjectId) %>% count()
trial <- subset(trial, trial$n > 35)
Data <- subset(Data,Data$subjectId %in% trial$subjectId)


# RT is positive if left food item chosen, negative if right food item chosen
idx = which(Data$choice==0)

Data$RT <- Data$rt
Data$RT[idx] = Data$rt[idx] * -1
Data$subject <- Data$subjectId

# scale value of health, taste and umami difference

Data <- Data %>% group_by(subjectId) %>%  mutate(hd = scale(healthrating)[,1])
Data <- Data %>% group_by(subjectId) %>%  mutate(td = scale(tasterating)[,1])
Data <- Data %>% group_by(subjectId) %>%  mutate(td = scale(umamirating)[,1])

hd = Data$healthrating
td = Data$tasterating
td = Data$umamirating

idxP = as.numeric(ordered(Data$subjectId)) #makes a sequentially numbered subj index
rtpos = Data$rt


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
sourceCpp("src/ddm_cpp_M4.cpp") 

# Parallel cross-validation
cat("Starting parallel cross-validation with", n_cores, "cores...\n")
fold_results <- foreach(k = 1:5, .combine = list, .multicombine = TRUE, .packages = c("runjags", "R2jags", "coda")) %dopar% {
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
                            time_umami.mu=0.1, time_umami.pr=0.05,
                            bias.mu=0.4,
                            bias.kappa=1, y_pred=y,  .RNG.name="base::Super-Duper", .RNG.seed=99999))


inits2 <- dump.format(list( alpha.mu=2.2, 
                            alpha.pr=0.05,  theta.mu=0.01,
                            theta.pr=0.05, b1.mu=0.3, b1.pr=0.05, b2.mu=0.1, b2.pr=0.05, b3.mu=0.1, b3.pr=0.05,
                            time_umami.mu=0.2, time_umami.pr=0.001,
                            bias.mu=0.4,
                            bias.kappa=1, y_pred=y,  .RNG.name="base::Wichmann-Hill", .RNG.seed=1234))

inits1 <- dump.format(list( alpha.mu=2.4,   
                            alpha.pr=0.05, theta.mu=0.15,
                            theta.pr=0.05, b1.mu=0.1, b1.pr=0.05, b2.mu=0.05, b2.pr=0.05, b3.mu=0.05, b3.pr=0.05,
                            time_umami.mu=0, time_umami.pr=0.01,
                            bias.mu=0.4,
                            bias.kappa=1, y_pred=y, .RNG.name="base::Mersenne-Twister", .RNG.seed=6666 ))
      
      monitor = c(
        "alpha.mu","theta.mu",
        "b1.mu","b2.mu","b3.mu", "bias.mu",
        "time_umami.mu",
        "b1.p","b2.p","b3.p",
        "time_umami.p",
        "theta.p", 
        "bias",
        "alpha.p")
      
      model = "src/ddm_priors_M4.txt"
      
     fit <- run.jags(model=model, 
                          monitor=monitor, data=dat, n.chains=3, inits=c(inits1,inits2, inits3), 
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
      b3_samples    <- get_param_matrix("b3")
      theta_samples <- get_param_matrix("theta")
      bias_samples  <- get_param_matrix_bias("bias")
      time_umami_samples <- get_param_matrix("time_umami")
      
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
          
          subframe <- data.frame(cbind(alpha_samples[,subj_idx],b1_samples[,subj_idx],b2_samples[,subj_idx],b3_samples[,subj_idx],theta_samples[,subj_idx],bias_samples[,subj_idx],time_umami_samples[,subj_idx]))
          colnames(subframe) <- c("alpha","b1","b2","b3","theta","bias","time_umami")
          hd <- subj_trials$hd[i]
          td <- subj_trials$td[i]
          pd <- subj_trials$pd[i]
          
          # Extract all posterior samples of parameters for this subject
          sim_samples <- numeric(n_iter)
          
          for (iter in 1:n_iter) {
            b1=subframe$b1[iter]
            b2=subframe$b2[iter]
            b3=subframe$b3[iter]
            alpha=subframe$alpha[iter]
            theta=subframe$theta[iter]
            bias=subframe$bias[iter]
            tIn_p=subframe$time_umami[iter]
            sim_samples[iter] <- ddm3_parallel(d_v = b1,d_h = b2,d_p = b3,thres = alpha,nDT = theta,bias = bias,vd =td ,hd =hd,pd =pd,tIn_p=tIn_p,N=1,sd_n=1)
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
        test_subject = subject_diffs
      )
      
      # Save individual fold result
      a <- data.frame(fold_result)
      a$subid <- rownames(a)
      write.csv(a, paste0("results/Sept_delayumami_cross_validated_fold_", k, ".csv"))
      
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

write.csv(s,"results/M4_cross_validated_folds.csv")
