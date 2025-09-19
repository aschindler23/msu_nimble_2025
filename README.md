# msu_nimble_2025
Introductary demo for running Bayesian models in NIMBLE with comparisons to JAGS

## Contents
- NIMBLE_demo_ppt.pdf: background info about NIMBLE with comparisons to JAGS
- `jags_model_parallel_3chains_ex.R`: example GLMM in JAGS
- `nimble_model_1chain_ex.R`: same GLMM example, but with one MCMC chain in NIMBLE
- `nimble_model_parallel_3chains_ex.R`: same GLMM example, but now running three MCMC chains in parallel
- `restart_nimble_run_ex.R`: same GLMM example but running one MCMC chain, saving progress, restarting R session, and continuing MCMC sampling from previous run
- `restart_nimble_run_parallel_ex.R`: same GLMM example but running three MCMC chains in parallel, saving progress, restarting R session, and continuing MCMC sampling from previous run
- `nimble_model_parallel_3chains_HMC_ex.R`: same GLMM example, but now running three MCMC chains in parallel using HMC NUTS
- `nimble_restart_functions.R`: functions from https://danielturek.github.io/public/saveMCMCstate/saveMCMCstate.html to save the model state and internal state of an MCMC algorithm to disk
