
### Poisson GLMM example in NIMBLE (three parallel chains) using HMC NUTS ###
# Author: A. Schindler

# load packages ----------------------------------------------------------------
library(nimble)
library(nimbleHMC)  # for running HMC NUTS
library(coda)       # for arranging results into MCMC objects
library(MCMCvis)    # for traceplots and summary statistics
library(parallel)   # for running parallel chains

# simulate data ----------------------------------------------------------------
set.seed(0)

# number data points
nobs <- 500

# number of sites
nsites <- 25

# site vector
site <- sample(1:nsites, nobs, replace = T)

# covariate data
x <- rnorm(nobs, 0, 1)

# simulate random site effect
epsilon <- rnorm(nsites, 0, 1)

# intercept
alpha <- 2

# slope
beta <- 0.3

# expected count
lambda <- rep(NA, nobs)
for(i in 1:nobs){
  lambda[i] <- exp(alpha + beta * x[i] + epsilon[site[i]])
}

# observed count
y <- rpois(nobs, lambda)

# inputs -----------------------------------------------------------------------
# write model code
model_code <- nimbleCode({
  ## priors
  ### fixed effects
  alpha ~ dnorm(0, sd = 10)
  beta ~ dnorm(0, sd = 10)
  
  ### random effects
  for(t in 1:nsites){
    epsilon[t] ~ dnorm(0, sd = sigma)
  }
  sigma ~ dunif(0, 10)
  
  for(i in 1:nobs){  
    ## likelihood
    log(lambda[i]) <- alpha + beta * x[i] + epsilon[site[i]]
    y[i] ~ dpois(lambda[i])
  }
})

# set constants
nimble_constants <- list(nobs = nobs, 
                         site = site,
                         nsites = nsites)

# set data
nimble_data <- list(y = y,
                    x = x)

# initial values
inits_function <- function(){list(
  alpha = rnorm(1, 0, 1),
  beta = rnorm(1, 0, 1),
  sigma = runif(1, 0, 1),
  epsilon = rnorm(nsites, 0, 1),
  lambda = y
)}

inits <- inits_function()

# create model object and check for errors -------------------------------------
model <- nimbleModel(code = model_code,
                     constants = nimble_constants,
                     data = nimble_data,
                     inits = inits,
                     calculate = F) # can set this to F when checking the model,
                                    # since we calculate again below

# check if fully initialized
model$initializeInfo()

# check for any errors
model$check()

# calculate logProb (should return a number)
model$calculate()

# run NIMBLE model in parallel -------------------------------------------------
# number of chains
nc <- 3

# set up cluster
cl <- makeCluster(nc)

# export inputs to the cluster
clusterExport(cl, c("model_code", "nimble_constants", "nimble_data"))

for(j in seq_along(cl)){
  nimble_inits <- inits_function() 
  clusterExport(cl[j], "nimble_inits")
}

# run MCMC in parallel with HMC NUTS
start_time <- Sys.time()
out <- clusterEvalQ(cl, {
  # load packages on cluster
  library(nimble)
  library(nimbleHMC)
  library(coda)
  library(parallel)
  
  # build model
  model <- nimbleModel(code = model_code, 
                       constants = nimble_constants,  
                       data =  nimble_data, 
                       inits = nimble_inits,
                       buildDerivs = T) # automatic derivatives needed
  
  # configure MCMC
  mcmc_Conf  <- configureHMC(model)
  
  # build MCMC
  modelMCMC  <- buildMCMC(mcmc_Conf)
  
  # compile model
  Cmodel     <- compileNimble(model)
  
  # compile MCMC
  CmodelMCMC <- compileNimble(modelMCMC, project = model)
  
  # run model
  CmodelMCMC$run(niter = 8000,
                 nburnin = 3000,
                 thin = 1)
  
  # extract samples
  samples <- as.matrix(CmodelMCMC$mvSamples)
  
  return(list(
    samples = samples
  ))
})
end_time <- Sys.time()
(run_time <- end_time - start_time)
stopCluster(cl)

# save results, traceplots, and results summary --------------------------------
# save output
file_heading <- "nimble_parallel_3chains_HMC_ex"
save(out, file = paste0(file_heading, "_results.RData"))

# convert to MCMC list
samples    <- list(chain1 = out[[1]][[1]], 
                   chain2 = out[[2]][[1]], 
                   chain3 = out[[3]][[1]])

mcmc_list <- as.mcmc.list(lapply(samples, mcmc))

# assess convergence
rhat <- gelman.diag(mcmc_list, multivariate = F)
rhat <- unlist(rhat$psrf)
rhat[which(rhat[,1] > 1.1), ]

# traceplots
MCMCtrace(mcmc_list, type = "trace",
          filename = paste0(file_heading, "_traceplots.pdf"))

# function to calculate percent of posterior above or below 0
prob_above_below_0 <- function(x){
  cdf <- ecdf(x)
  prob_below_0 <- cdf(0)
  prob_above_0 <- 1 - prob_below_0
  if(median(x) < 0){
    return(prob_below_0)
  } else{
    return(prob_above_0)
  }
}

# calculate summary statistics
sum_stats <- MCMCsummary(mcmc_list, func = function(x) prob_above_below_0(x), func_name = "f")

# save summary statistics
write.csv(sum_stats, paste0(file_heading, "_summary.csv"))

