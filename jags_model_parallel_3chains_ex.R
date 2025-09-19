
### Poisson GLMM example in JAGS (three parallel chains) ###
# Author: A. Schindler

# load packages ----------------------------------------------------------------
library(jagsUI)
library(coda)

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

# run model in JAGS ------------------------------------------------------------
# model code
sink("jags_model_ex.txt")
cat("
model {
  ## priors
  ### fixed effects
  alpha ~ dnorm(0, 0.01)
  beta ~ dnorm(0, 0.01)
  
  ### random effects
  for(t in 1:nsites){
    epsilon[t] ~ dnorm(0, tau)
  }
  tau <- pow(sigma, -2)
  sigma ~ dunif(0, 10)
  
  for(i in 1:nobs){  
    ## likelihood
    log(lambda[i]) <- alpha + beta * x[i] + epsilon[site[i]]
    y[i] ~ dpois(lambda[i])
  }
}", fill = T)
sink()

# set data
jags_data <- list(
  nobs = nobs, 
  site = site,
  nsites = nsites,
  y = y,
  x = x
)

# set initial values
inits_function <- function(){list(
  alpha = rnorm(1, 0, 1),
  beta = rnorm(1, 0, 1),
  sigma = runif(1, 0, 1)
)}

# parameters to monitor
params <- c("alpha",
            "beta",
            "sigma")

# run model
start_time <- Sys.time()
out <- jags(jags_data, inits_function, params, "jags_model_ex.txt",
            n.thin = 10, n.chains = 3, n.burnin = 30000, n.iter = 80000,
            parallel = T)
end_time <- Sys.time()
(run_time <- end_time - start_time)

# save results, traceplots, and results summary --------------------------------
# save output
file_heading <- "jags_parallel_3chains_ex"
save(out, file = paste0(file_heading, "_results.RData"))

# assess convergence
range(out$Rhat, na.rm = T)

pdf(paste0(file_heading, "_traceplots.pdf"))
jagsUI::traceplot(out)
dev.off()

# save summary statistics
write.csv(out$summary, paste0(file_heading, "_summary.csv"))

