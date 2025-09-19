
### Functions to "save" the model state and internal state of an MCMC algorithm to disk
### from: https://danielturek.github.io/public/saveMCMCstate/saveMCMCstate.html

getStateVariableNames <- function(samplerDef) {
  resetMethod <- body(samplerDef$reset)
  stateVars <- character()
  if(resetMethod[[1]] != '{') stop('something wrong')
  numLines <- length(resetMethod)
  for(i in 1:numLines) {
    if(i == 1) next
    thisLine <- resetMethod[[i]]
    if(thisLine[[1]] == '<<-') {
      LHS <- thisLine[[2]]
      if(!is.name(LHS)) stop('haven\'t dealt with non-name-LHS case yet')
      stateVars <- c(stateVars, as.character(LHS))
    }
    if('my_calcAdaptationFactor' %in% all.names(thisLine)) {
      stateVars <- c(stateVars, 'my_calcAdaptationFactor')
    }
  }
  setupMethod <- body(samplerDef$setup)
  if('empirSamp' %in% all.names(setupMethod)) stateVars <- c(stateVars, 'empirSamp')
  return(stateVars)
}

getModelState <- function(model) {
  modelVarNames <- model$getVarNames()
  modelVarValuesList <- vector('list', length(modelVarNames))
  names(modelVarValuesList) <- modelVarNames
  for(var in modelVarNames) {
    modelVarValuesList[[var]] <- model[[var]]
  }
  return(modelVarValuesList)
}

getMCMCstate <- function(conf, mcmc) {
  stateVarNamesList <- vector('list', length(conf$samplerConfs))
  mcmcStateValuesList <- vector('list', length(conf$samplerConfs))
  for(i in seq_along(conf$samplerConfs)) {
    samplerDef <- conf$getSamplerDefinition(i)
    theseStateNames <- getStateVariableNames(samplerDef)
    theseStateValuesList <- vector('list', length(theseStateNames))
    names(theseStateValuesList) <- theseStateNames
    for(j in seq_along(theseStateNames)) {
      if(is.nf(mcmc)) {
        if(theseStateNames[j] == 'my_calcAdaptationFactor') {
          theseStateValuesList[[j]] <- list(timesAdapted = mcmc$samplerFunctions$contentsList[[i]][[theseStateNames[j]]]$timesAdapted,
                                            gamma1 = mcmc$samplerFunctions$contentsList[[i]][[theseStateNames[j]]]$gamma1)
        } else
          theseStateValuesList[[j]] <- mcmc$samplerFunctions$contentsList[[i]][[theseStateNames[j]]]
      }
      if(is.Cnf(mcmc)) {
        if(theseStateNames[j] == 'my_calcAdaptationFactor') {
          theseStateValuesList[[j]] <- list(timesAdapted = valueInCompiledNimbleFunction(mcmc$Robject$samplerFunctions$contentsList[[i]]$my_calcAdaptationFactor$.CobjectInterface, 'timesAdapted'),
                                            gamma1 = valueInCompiledNimbleFunction(mcmc$Robject$samplerFunctions$contentsList[[i]]$my_calcAdaptationFactor$.CobjectInterface, 'gamma1'))
        } else
          theseStateValuesList[[j]] <- valueInCompiledNimbleFunction(mcmc$samplerFunctions[[i]], theseStateNames[j])
      }
    }
    mcmcStateValuesList[[i]] <- theseStateValuesList
  }
  return(mcmcStateValuesList)
}

getWAICstate <- function(mcmc) {
  waicStateNames <- c('delta1pWAICmat', 'delta2pWAICmat', 'finalized', 'logProbMat', 'lppdCurSumMat', 'lppdSumMaxMat',
                      'mcmcIter', 'meanpWAICmat', 'sspWAICmat')
  waicStateList <- vector('list', length(waicStateNames))
  names(waicStateList) <- waicStateNames
  for(nm in waicStateNames) {
    if(is.nf(mcmc)) {
      waicStateList[[nm]] <- mcmc$waicFun[[1]][[nm]]
    }
    if(is.Cnf(mcmc)) {
      waicStateList[[nm]] <- valueInCompiledNimbleFunction(mcmc$waicFun[[1]], nm)
    }
  }
  return(waicStateList)
}

setModelState <- function(model, modelState) {
  modelVarNames <- model$getVarNames()
  if(!identical(sort(modelVarNames), sort(names(modelState)))) stop('saved model variables don\'t agree')
  for(var in modelVarNames) {
    model[[var]] <- modelState[[var]]
  }
  invisible(model$calculate())
}

setMCMCstate <- function(conf, mcmc, mcmcState) {
  if(length(mcmcState) != length(conf$samplerConfs)) stop('saved mcmc samplers don\'t agree')
  for(i in seq_along(conf$samplerConfs)) {
    theseStateValuesList <- mcmcState[[i]]
    for(j in seq_along(theseStateValuesList)) {
      samplerStateName <- names(theseStateValuesList)[j]
      if(is.nf(mcmc)) {
        if(samplerStateName == 'my_calcAdaptationFactor') {
          mcmc$samplerFunctions$contentsList[[i]]$my_calcAdaptationFactor$timesAdapted <- theseStateValuesList[[j]]$timesAdapted
          mcmc$samplerFunctions$contentsList[[i]]$my_calcAdaptationFactor$gamma1 <- theseStateValuesList[[j]]$gamma1
        } else {
          mcmc$samplerFunctions$contentsList[[i]][[samplerStateName]] <- theseStateValuesList[[samplerStateName]]
        }
      }
      if(is.Cnf(mcmc)) {
        if(samplerStateName == 'my_calcAdaptationFactor') {
          invisible(valueInCompiledNimbleFunction(mcmc$Robject$samplerFunctions$contentsList[[i]]$my_calcAdaptationFactor$.CobjectInterface, 'timesAdapted', theseStateValuesList[[j]]$timesAdapted))
          invisible(valueInCompiledNimbleFunction(mcmc$Robject$samplerFunctions$contentsList[[i]]$my_calcAdaptationFactor$.CobjectInterface, 'gamma1', theseStateValuesList[[j]]$gamma1))
        } else {
          invisible(valueInCompiledNimbleFunction(mcmc$samplerFunctions[[i]], samplerStateName, theseStateValuesList[[samplerStateName]]))
        }
      }
    }
  }
}

setWAICstate <- function(mcmc, waicState) {
  for(nm in names(waicState)) {
    if(is.nf(mcmc)) {
      mcmc$waicFun[[1]][[nm]] <- waicState[[nm]]
    }                       
    if(is.Cnf(mcmc)) {
      invisible(valueInCompiledNimbleFunction(mcmc$waicFun[[1]], nm, waicState[[nm]]))
    }
  }
}