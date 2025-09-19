
### example of how to restart multiple parallel chains of a NIMBLE MCMC algorithm 
#    from where it left off, assuming you have to end your R session between runs
# adapted from: https://danielturek.github.io/public/restartingMCMC/restartingMCMC.html 
# simplified with functions from: https://danielturek.github.io/public/saveMCMCstate/saveMCMCstate.html

### Poisson GLMM example in NIMBLE (one chain) ###
# Author: A. Schindler

# load packages ----------------------------------------------------------------
library(nimble)
library(parallel)
source("nimble_restart_functions.R")

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
                     calculate = F)

# check if fully initialized
model$initializeInfo()

# check for any errors
model$check()

# calculate logProb (should return a number)
model$calculate()

# run NIMBLE model in parallel -------------------------------------------------
# select number of chains
nc <- 3

# set seed for model run (FOR ILLUSTRATIVE PURPOSES ONLY, REMOVE FOR REAL MODELS)
seed <- c(1, 2, 3)

# set up cluster
cl <- makeCluster(nc)

clusterExport(cl, c("model_code", "nimble_data", "nimble_constants", "inits", 
                    "getMCMCstate", "getModelState", "getStateVariableNames",
                    "getWAICstate", "setMCMCstate", "setModelState", "setWAICstate"))

for(j in seq_along(cl)){
  nimble_inits <- inits_function() 
  nimble_seed <- seed[j] # FOR ILLUSTRATIVE PURPOSES ONLY, REMOVE FOR REAL MODELS
  clusterExport(cl[j], c("nimble_inits", "nimble_seed"))
}

# run MCMC in parallel
start_time <- Sys.time()
out <- clusterEvalQ(cl, {
  library(nimble)
  library(coda)
  library(parallel)

  # build model
  model <- nimbleModel(code = model_code, 
                       constants = nimble_constants,  
                       data =  nimble_data, 
                       inits = nimble_inits)

  # create and configure MCMC algorithm
  useWAIC <- F # can turn on if you want nimble to calculate WAIC
  mcmc_Conf <- configureMCMC(model, enableWAIC = useWAIC)
  
  # build the MCMC
  modelMCMC <- buildMCMC(mcmc_Conf)
  
  # compile model and MCMC -------------------------------------------------------
  # compile model
  Cmodel <- compileNimble(model)
  
  # compile MCMC
  CmodelMCMC <- compileNimble(modelMCMC, project = model)
  
  # set seed for model run (FOR ILLUSTRATIVE PURPOSES ONLY, REMOVE FOR REAL MODELS)
  set.seed(nimble_seed)

  # run model ------------------------------------------------------------------
  CmodelMCMC$run(niter = 80000,
                 nburnin = 30000,
                 thin = 10)

  ### What to save
  ## 1. MCMC samples
  ## 2. values of all non-data stochastic nodes in the model (model state)
  ## 3. internal state variables from all the MCMC samplers
  ## 4. R's random number seed
  ## 5. internal state variables for WAIC calculation (if used)

  # 1. extract samples
  samples <- as.matrix(CmodelMCMC$mvSamples)

  # 2. extract model state 
  modelState <- getModelState(Cmodel)
  
  # 3. extract internal state variables
  mcmcState <- getMCMCstate(mcmc_Conf, CmodelMCMC)

  # 4. extract R's random number seed
  seed <- .Random.seed
  
  # 5. extract WAIC state if used
  if(useWAIC){
    waicState <- getWAICstate(CmodelMCMC)
  } else{
    waicState <- NULL
  }

  return(list(
    samples = samples,
    modelState = modelState,
    mcmcState = mcmcState,
    seed = seed,
    waicState = waicState
  ))
})
end_time <- Sys.time()
(run_time <- end_time - start_time)
stopCluster(cl)

# save state files
fileName <- "nimble_restart_files.RData"
save(out, file = fileName)

### Restart R to demonstrate a new R session
# load packages ----------------------------------------------------------------
library(nimble)
library(parallel)
source("nimble_restart_functions.R")

# load results
fileName <- "nimble_restart_files.RData"
load(fileName)

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
                     calculate = F)

# check if fully initialized
model$initializeInfo()

# check for any errors
model$check()

# calculate logProb (should return a number)
model$calculate()

# resume model run -------------------------------------------------------------
# select number of chains
nc <- 3

# set up cluster
cl <- makeCluster(nc)

clusterExport(cl, c("model_code", "nimble_data", "nimble_constants", "inits", 
                    "getMCMCstate", "getModelState", "getStateVariableNames",
                    "getWAICstate", "setMCMCstate", "setModelState", "setWAICstate"))

for(j in seq_along(cl)){
  nimble_inits <- inits_function() 
  modelState <- out[[j]]$modelState
  mcmcState <- out[[j]]$mcmcState
  seed <- out[[j]]$seed
  waicState <- out[[j]]$waicState
  clusterExport(cl[j], c("nimble_inits", "modelState", "mcmcState", "seed", "waicState"))
}

# run MCMC in parallel
start_time <- Sys.time()
out_update <- clusterEvalQ(cl, {
  library(nimble)
  library(coda)
  library(parallel)
  
  # build model
  model <- nimbleModel(code = model_code, 
                       constants = nimble_constants,  
                       data =  nimble_data, 
                       inits = nimble_inits)
  
  # create and configure MCMC algorithm
  useWAIC <- F # can turn on if you want nimble to calculate WAIC
  mcmc_Conf <- configureMCMC(model, enableWAIC = useWAIC)
  
  # build the MCMC
  modelMCMC <- buildMCMC(mcmc_Conf)
  
  # compile model and MCMC -------------------------------------------------------
  # compile model
  Cmodel <- compileNimble(model)
  
  # compile MCMC
  CmodelMCMC <- compileNimble(modelMCMC, project = model)

  # restore previous state
  setModelState(Cmodel, modelState)
  setMCMCstate(mcmc_Conf, CmodelMCMC, mcmcState)
  if(useWAIC){
    setWAICstate(CmodelMCMC, waicState)
  }
  .Random.seed <- seed  

  # resume MCMC run
  CmodelMCMC$run(niter = 50000,
            thin = 10,
            reset = FALSE, 
            resetWAIC = FALSE)

  # extract samples
  samples <- as.matrix(CmodelMCMC$mvSamples)
  
  # extract model state and internal state of MCMC 
  modelState <- getModelState(Cmodel)
  mcmcState <- getMCMCstate(mcmc_Conf, CmodelMCMC)
  
  # extract R's random number seed
  seed <- .Random.seed
  
  # extract WAIC state if used
  if(useWAIC){
    waicState <- getWAICstate(Cmcmc)
  } else{
    waicState <- NULL
  }
  
  return(list(
    samples = samples,
    modelState = modelState,
    mcmcState = mcmcState,
    seed = seed,
    waicState = waicState
  ))
})
end_time <- Sys.time()
(run_time <- end_time - start_time)
stopCluster(cl)

# extract original samples
samples <- list(chain1 = out[[1]]$samples, 
                chain2 = out[[2]]$samples, 
                chain3 = out[[3]]$samples)

# extract new samples
samples_continued <- list(chain1 = out_update[[1]]$samples[-1,], 
                          chain2 = out_update[[2]]$samples[-1,], 
                          chain3 = out_update[[3]]$samples[-1,])

# look at samples from first run, and from second
chain1_all <- rbind(samples$chain1, samples_continued$chain1)
chain2_all <- rbind(samples$chain2, samples_continued$chain2)
chain3_all <- rbind(samples$chain3, samples_continued$chain3)

chain1_all[4995:5005,]
chain2_all[4995:5005,]
chain3_all[4995:5005,]

### Verify this is the same as one long run
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
                     calculate = F)

# check if fully initialized
model$initializeInfo()

# check for any errors
model$check()

# calculate logProb (should return a number)
model$calculate()

# run NIMBLE model in parallel -------------------------------------------------
# select number of chains
nc <- 3

# set seed for model run (FOR ILLUSTRATIVE PURPOSES ONLY, REMOVE FOR REAL MODELS)
seed <- c(1, 2, 3)

# set up cluster
cl <- makeCluster(nc)

clusterExport(cl, c("model_code", "nimble_data", "nimble_constants", "inits", 
                    "getMCMCstate", "getModelState", "getStateVariableNames",
                    "getWAICstate", "setMCMCstate", "setModelState", "setWAICstate"))

for(j in seq_along(cl)){
  nimble_inits <- inits_function() 
  nimble_seed <- seed[j] # FOR ILLUSTRATIVE PURPOSES ONLY, REMOVE FOR REAL MODELS
  clusterExport(cl[j], c("nimble_inits", "nimble_seed"))
}

# run MCMC in parallel
start_time <- Sys.time()
out_full <- clusterEvalQ(cl, {
  library(nimble)
  library(coda)
  library(parallel)
  
  # build model
  model <- nimbleModel(code = model_code, 
                       constants = nimble_constants,  
                       data =  nimble_data, 
                       inits = nimble_inits)
  
  # create and configure MCMC algorithm
  useWAIC <- F # can turn on if you want nimble to calculate WAIC
  mcmc_Conf <- configureMCMC(model, enableWAIC = useWAIC)
  
  # build the MCMC
  modelMCMC <- buildMCMC(mcmc_Conf)
  
  # compile model and MCMC -------------------------------------------------------
  # compile model
  Cmodel <- compileNimble(model)
  
  # compile MCMC
  CmodelMCMC <- compileNimble(modelMCMC, project = model)
  
  # set seed for model run (FOR ILLUSTRATIVE PURPOSES ONLY, REMOVE FOR REAL MODELS)
  set.seed(nimble_seed)
  
  # run model ------------------------------------------------------------------
  CmodelMCMC$run(niter = 130000,
                 nburnin = 30000,
                 thin = 10)
  
  ### What to save
  ## 1. MCMC samples
  ## 2. values of all non-data stochastic nodes in the model (model state)
  ## 3. internal state variables from all the MCMC samplers
  ## 4. R's random number seed
  ## 5. internal state variables for WAIC calculation (if used)
  
  # 1. extract samples
  samples <- as.matrix(CmodelMCMC$mvSamples)
  
  # 2. extract model state 
  modelState <- getModelState(Cmodel)
  
  # 3. extract internal state variables
  mcmcState <- getMCMCstate(mcmc_Conf, CmodelMCMC)
  
  # 4. extract R's random number seed
  seed <- .Random.seed
  
  # 5. extract WAIC state if used
  if(useWAIC){
    waicState <- getWAICstate(CmodelMCMC)
  } else{
    waicState <- NULL
  }
  
  return(list(
    samples = samples,
    modelState = modelState,
    mcmcState = mcmcState,
    seed = seed,
    waicState = waicState
  ))
})
end_time <- Sys.time()
(run_time <- end_time - start_time)
stopCluster(cl)


# extract samples
samples_full <- list(chain1 = out_full[[1]]$samples, 
                     chain2 = out_full[[2]]$samples, 
                     chain3 = out_full[[3]]$samples)

# look at samples from first run, and from second
chain1_all[4995:5005,]
samples_full$chain1[4995:5005,]

chain2_all[4995:5005,]
samples_full$chain2[4995:5005,]

chain3_all[4995:5005,]
samples_full$chain3[4995:5005,]

