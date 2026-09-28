# Adapted from code originally written by Blair Shevlin for "Negative affect
# influences the computations underlying food choice in bulimia nervosa".

# This version fits Time-Varying DDM with Taste and Health only (thonly)

# Packages required
required_packages <- c("DEoptim", 
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
  "coda")

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

#master$fat.e <- as.factor(master$fat.e)

master$choice.bin <- ifelse(master$choice.rating > 5,1,0)
master$choice.bin[master$choice.rating==5] <- NA
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


#data$tasterating <- data$vt_l - data$vt_r
#data$healthrating <- data$vh_l - data$vh_r

#do not need to standardize within individual OR compare with neutral item for difference since
# the neutral item is constant
data$tasterating <- data$vt_l - data$vt_r
data$healthrating <- data$vh_l - data$vh_r
data$umamirating <- data$Umami_og - data$ref_p

fi <- data
name <- "mTurk"
ALL <- fi %>% dplyr::select(subjectId, tasterating, healthrating, umamirating,rt, choice)
ALL[,2:6] <- sapply(ALL[,2:6], as.numeric)

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
Data <- Data %>% group_by(subjectId) %>%  mutate(pd = scale(umamirating)[,1])

hd = Data$hd
td = Data$td
pd = Data$pd

idxP = as.numeric(ordered(Data$subjectId)) #makes a sequentially numbered subj index
rtpos = Data$rt

subs = Data %>% dplyr::select(subjectId) %>% unique()
write.csv(subs,"results/new_Sept_TV_thonly_sublist.csv")

y= Data$RT
N = length(y)
ns = length(unique(idxP))        


set.seed(42)


# Define number of folds (or train/test split)
K <- 5
Data$fold <- rep(1:K, length.out = nrow(Data))

# Simulation code (load once before parallel processing)
sourceCpp("src/ddm_cpp_M3.cpp")

# Parallel cross-validation
cat("Fitting on FULL data (no CV)...\n")
fold_results <- vector("list", 5)
  train_data <- Data   # FULL data: no fold split
  idxP = as.numeric(ordered(Data$subjectId)) 
  y=Data$RT
  N=length(Data$RT)
  idxP=idxP
  hd=Data$hd
  td=Data$td
  rtpos = Data$rt
  ns = length(unique(idxP))
  
  
  jags_data <- dump.format(list(N=N, y=y, idxP=idxP, hd=hd, td=td, rt=rtpos, ns=ns))
  
      
      inits3 <- dump.format(list( alpha.mu=2, time.mu=0.1, 
                       alpha.pr=0.5, time.pr= 0.5, theta.mu=0.1,
                       theta.pr=0.05,  b1.mu=0.3, b1.pr=0.05, b2.mu=0.01, b2.pr=0.05, 
                       bias.mu=0.4,
                       bias.kappa=1, y_pred=y,  .RNG.name="base::Super-Duper", .RNG.seed=99999))
      
      
      inits2 <- dump.format(list( alpha.mu=2.2, time.mu=-0.1, 
                       alpha.pr=0.05, time.pr= 0.05, theta.mu=0.01,
                       theta.pr=0.05, b1.mu=0.3, b1.pr=0.05, b2.mu=0.1, b2.pr=0.05, 
                       bias.mu=0.4,
                       bias.kappa=1, y_pred=y,  .RNG.name="base::Wichmann-Hill", .RNG.seed=1234))
      
      inits1 <- dump.format(list( alpha.mu=2.4, time.mu=0,  
                       alpha.pr=0.05, time.pr= 0.05, theta.mu=0.15,
                       theta.pr=0.05, b1.mu=0.1, b1.pr=0.05, b2.mu=0.05, b2.pr=0.05, 
                       bias.mu=0.4,
                       bias.kappa=1, y_pred=y, .RNG.name="base::Mersenne-Twister", .RNG.seed=6666 ))
      
      monitor = c(
        "time.mu","alpha.mu","theta.mu",
        "b1.mu","b2.mu","bias.mu",
        "time.p",
        "b1.p","b2.p", 
        "theta.p", 
        "bias",
        "alpha.p","e.p.tau")
      
      model = "src/ddm_priors_M3.txt"
      
      fit <- runjags::run.jags(model=model, 
                                   monitor=monitor, data=jags_data, n.chains=3, inits=c(inits1,inits2, inits3), 
                                   plots = FALSE, method="parallel", module="wiener", burnin=50000, sample=10000, thin=1)

# ---- save a thinned posterior (full object is ~1GB; 3000 draws is ample for WAIC/LOO + recovery) ----
samples <- as.mcmc.list(fit)
comb <- do.call(rbind, lapply(samples, as.matrix))
keep <- round(seq(1, nrow(comb), length.out = min(3000, nrow(comb))))
post <- comb[keep, , drop = FALSE]

dir.create("results/fits", showWarnings = FALSE, recursive = TRUE)
saveRDS(list(model = "M3", post = post,
             n_obs = nrow(Data), subjects = levels(droplevels(factor(Data$subjectId))),
             Data = as.data.frame(Data)),
        file = "results/fits/fullfit_M3.rds")
write.csv(as.data.frame(Data), file = "results/fits/Data_M3.csv", row.names = FALSE)
cat("\n=== M3 full fit saved: results/fits/fullfit_M3.rds ===\n")
cat("posterior draws kept:", nrow(post), " params:", ncol(post), "\n")
