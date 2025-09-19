
### example of how to restart a NIMBLE MCMC algorithm from where it left off, 
#    assuming you have to end your R session between runs
# adapted from: https://danielturek.github.io/public/restartingMCMC/restartingMCMC.html 
# simplified with functions from: https://danielturek.github.io/public/saveMCMCstate/saveMCMCstate.html

### Poisson GLMM example in NIMBLE (one chain) ###
# Author: A. Schindler

# load packages ----------------------------------------------------------------
library(nimble)
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
                     inits = inits)

# check if fully initialized
model$initializeInfo()

# check for any errors
model$check()

# calculate logProb (should return a number)
model$calculate()

# build MCMC -------------------------------------------------------------------
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

# run model (one chain) --------------------------------------------------------
set.seed(0) # for demonstration only, don't want to typically do this
CmodelMCMC$run(niter = 80000,
               nburnin = 30000,
               thin = 10)

# extract samples --------------------------------------------------------------

### What to save
## 1. MCMC samples
## 2. values of all non-data stochastic nodes in the model (model state)
## 3. internal state variables from all the MCMC samplers
## 4. R's random number seed
## 5. internal state variables for WAIC calculation (if used)

#  1. extract samples
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
}

# save state files
fileName <- "nimble_restart_files.RData"

if(useWAIC){
  save(samples, modelState, mcmcState, waicState, seed, file = fileName)
} else{
  save(samples, modelState, mcmcState, seed, file = fileName)
}

### Start new R session ###
# load packages ----------------------------------------------------------------
library(nimble)
source("nimble_restart_functions.R")

# load results
load("nimble_restart_files.RData")

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
                     inits = inits)

# check if fully initialized
model$initializeInfo()

# check for any errors
model$check()

# calculate logProb (should return a number)
model$calculate()

# build MCMC -------------------------------------------------------------------
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
               reset = F,
               resetWAIC = F)

# extract new samples
samples_continued <- as.matrix(CmodelMCMC$mvSamples)[-1,]

# look at samples from first run, and from second
samples_all <- rbind(samples, samples_continued)
dim(samples_all)

samples_all[4995:5005,]

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
                     inits = inits)

# check if fully initialized
model$initializeInfo()

# check for any errors
model$check()

# calculate logProb (should return a number)
model$calculate()

# build MCMC -------------------------------------------------------------------
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

# run model (one chain) --------------------------------------------------------
set.seed(0) # for demonstration only, don't want to typically do this
CmodelMCMC$run(niter = 130000,
               nburnin = 30000,
               thin = 10)

# extract samples
samples_all_long <- as.matrix(CmodelMCMC$mvSamples)
dim(samples_all_long)

# compare samples
samples_all[4995:5005,]
samples_all_long[4995:5005,]

