
### Poisson GLMM example in NIMBLE (one chain) ###
# Author: A. Schindler

# load packages ----------------------------------------------------------------
library(nimble)

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

# step 1: inputs ---------------------------------------------------------------
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

# step 2: create model object and check for errors -----------------------------
# create model (note, can set calculate = F to speed up making model object here for complex models)
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

# step 3: build MCMC -----------------------------------------------------------
# create and configure MCMC algorithm
mcmc_Conf <- configureMCMC(model)

# look at samplers if desired
mcmc_Conf$printSamplers()

# an example of how to change samplers - replace RW sampler with slice sampler
mcmc_Conf$removeSampler("sigma")
mcmc_Conf$addSampler("sigma", "slice")
mcmc_Conf$printSamplers()

# note: only parameters with priors are saved by default, but you can tell NIMBLE 
#       to monitor additional parameters if desired 
mcmc_Conf$printMonitors()
mcmc_Conf$addMonitors("lambda")

# build the MCMC
modelMCMC <- buildMCMC(mcmc_Conf)

# step 4: compile model and MCMC -----------------------------------------------
# compile model
Cmodel <- compileNimble(model)

# compile MCMC
CmodelMCMC <- compileNimble(modelMCMC, project = model)

# step 5: run model (one chain) ------------------------------------------------
CmodelMCMC$run(niter = 80000,
               nburnin = 30000,
               thin = 10)

# step 6: extract samples ------------------------------------------------------
samples <- as.matrix(CmodelMCMC$mvSamples)

# save samples
save(samples, file = "nimble_1chain_ex_samples.RData")

