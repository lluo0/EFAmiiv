

EFAmiiv_multi_s1 <- function(model = NULL,
                          data,
                          sigLevel = .05,
                          scalingCrit = "sargan+factorloading_R2",
                          correlatedErrors = NULL){

  ##get the initial model with variables on each factor saved separately
  modelparts <- strsplit(model, split = '\n', fixed = T)[[1]]
  num_factor <- length(modelparts)
  model_sep <- list()
  for(n in 1:num_factor){
    model_sep[[n]] <- strsplit(strsplit(modelparts, split = '=~', fixed = T)[[n]][2], split = '+', fixed = T)[[1]]
  }
  model_sep <- lapply(model_sep, function(x) gsub(" ", "",x))

  if(!is.null(correlatedErrors)){
    model <- paste0(model, '\n', correlatedErrors)
  }
  fit_initial <- miive(model = model, data = data, var.cov = T)

  # ##get bad variables
  # badvar <- getbadvar(fit_initial)


  ##get bad variables
  badvar <- getbadvar(fit_initial, sigLevel = sigLevel)

  ##crossload bad variables
  crossloadmodel_sep <- lapply(model_sep, function(x) unique(c(x, badvar)))

  ##crossload model
  crossloadmodel <- list()
  for(n in 1:num_factor){
    crossloadmodel[[n]] <- paste0('f', n, '=~', paste0(crossloadmodel_sep[[n]], collapse = '+'))
  }
  ##crossload fit
  crossloadfit <- miive(model = paste0(paste0(crossloadmodel, collapse = '\n'), '\n', correlatedErrors),
                        data, var.cov = T)

  ##detect unnecessary cross loadings
  badvar_crossload <- getbadvar_crossload(crossloadfit, num_fac = num_factor, badvar = badvar)

  ##remove unncessary cross loadings
  crossloadmodel_sepclean <- list()
  for(n in 1:num_factor){
    crossloadmodel_sepclean[[n]] <- setdiff(crossloadmodel_sep[[n]], badvar_crossload[[n]])
  }

  ##new crossload model
  crossloadmodel_clean <- list()
  for(n in 1:num_factor){
    crossloadmodel_clean[[n]] <- paste0('f', n, '=~', paste0(crossloadmodel_sepclean[[n]], collapse = '+'))
  }

  ##check if any variables not loaded on any factors!
  anyleftout <- setdiff(colnames(data),unique(unlist(crossloadmodel_sepclean)))
  ##load left out variable on the last factor if only one left out variable, otherwise need to create a new factor
  if(length(anyleftout) == 1){
    crossloadmodel_clean[[num_factor]] <- paste0(crossloadmodel_clean[[num_factor]], '+', anyleftout)
  }
  if(length(anyleftout) == 0){
    crossloadmodel_clean <- crossloadmodel_clean
  }
  if(length(anyleftout)>1){
    #setup stepPrev object to select scaling indiactor for this new factor
    stepPrev <- list()
    stepPrev$goodmodelpart <- crossloadmodel_clean
    stepPrev$badvar <- anyleftout
    stepPrev$num_factor <- num_factor
    stepPrev$correlatedErrors <- correlatedErrors
    scalingind <- select_scalingind_stepN(data, sigLevel, scalingCrit,
                                          stepPrev = stepPrev)
    crossloadmodel_clean[[num_factor+1]] <- paste0('f', num_factor+1, '=~',
                                                   scalingind, '+',
                                                   paste0(anyleftout[-which(anyleftout==scalingind)], collapse = '+'))
  }

  ##new fit
  newcrossloadfit <- miive(model = paste0(paste0(crossloadmodel_clean, collapse = '\n'), '\n', correlatedErrors), data, var.cov = T)

  ##get badvar
  newbadvar <- getbadvar(newcrossloadfit, sigLevel = sigLevel)


  ##decide if need next step
  nextstep <- 'no'
  if(length(newbadvar) > 1 | length(anyleftout) > 1){
    nextstep <- 'yes'
  }

  # if(!is.null(correlatedErrors)){
  #   finalmodel <- paste0(paste0(crossloadmodel_clean, collapse = '\n'), '\n', correlatedErrors)
  # }else{finalmodel <- paste0(crossloadmodel_clean, collapse = '\n')}

  ##create the final obj to return
  finalobj <- list(#model = finalmodel,
                   model = paste0(crossloadmodel_clean, collapse ='\n'),
                   fit = newcrossloadfit,
                   num_factor = num_factor,
                   # badvar = c(newbadvar, anyleftout), #need to combine the badvars and left out variables
                   badvar = newbadvar,
                   nextstep = nextstep,
                   correlatedErrors = correlatedErrors,
                   model_sep = crossloadmodel_sepclean,
                   goodmodelpart = crossloadmodel_clean)
  return(finalobj)
}

EFAmiiv_multi_s2 <- function(stepPrev,
                             data,
                             sigLevel =  .05,
                             scalingCrit = "sargan+factorloading_R2"){

  # ##create a new factor for bad variables
  # scalingindicator <- select_scalingind_stepN(data, scalingCrit = scalingCrit, stepPrev = stepPrev)
  # correlatedErrors <- stepPrev$correlatedErrors
  # ##create the new model
  # num_factor <- stepPrev$num_factor + 1
  # model_sep <- stepPrev$model_sep
  # model_sep[[num_factor]] <- c(scalingindicator, setdiff(stepPrev$badvar, scalingindicator))
  # model <- list()
  # for(n in 1:num_factor){
  #   model[[n]] <- paste0('f', n, '=~', paste0(model_sep[[n]], collapse = '+'))
  # }
  # ##get fit
  # fit <- miive(paste0(paste0(model, collapse = '\n'), '\n', correlatedErrors), data, var.cov = T)
  # ##get bad variables
  # badvar <- getbadvar(fit, sigLevel = sigLevel)
  model <- stepPrev$model
  fit <- stepPrev$fit
  badvar <- stepPrev$badvar
  correlatedErrors <- stepPrev$correlatedErrors

  ##get the initial model with variables on each factor saved separately
  modelparts <- strsplit(model, split = '\n', fixed = T)[[1]]
  num_factor <- length(modelparts)
  model_sep <- list()
  for(n in 1:num_factor){
    model_sep[[n]] <- strsplit(strsplit(modelparts, split = '=~', fixed = T)[[n]][2], split = '+', fixed = T)[[1]]
  }
  model_sep <- lapply(model_sep, function(x) gsub(" ", "",x))
  ##crossload bad variables
  crossloadmodel_sep <- lapply(model_sep, function(x) unique(c(x, badvar)))

  ##crossload model
  crossloadmodel <- list()
  for(n in 1:num_factor){
    crossloadmodel[[n]] <- paste0('f', n, '=~', paste0(crossloadmodel_sep[[n]], collapse = '+'))
  }
  ##crossload fit
  crossloadfit <- miive(model = paste0(paste0(crossloadmodel, collapse = '\n'), '\n', correlatedErrors),
                        data, var.cov = T)

  ##detect unnecessary cross loadings
  badvar_crossload <- getbadvar_crossload(crossloadfit, num_fac = num_factor, badvar = badvar)

  ##remove unncessary cross loadings
  crossloadmodel_sepclean <- list()
  for(n in 1:num_factor){
    crossloadmodel_sepclean[[n]] <- setdiff(crossloadmodel_sep[[n]], badvar_crossload[[n]])
  }

  ##new crossload model
  crossloadmodel_clean <- list()
  for(n in 1:num_factor){
    crossloadmodel_clean[[n]] <- paste0('f', n, '=~', paste0(crossloadmodel_sepclean[[n]], collapse = '+'))
  }

  ##check if any variables not loaded on any factors!
  anyleftout <- setdiff(colnames(data),unique(unlist(crossloadmodel_sepclean)))
  ##load left out variable on the last factor if only one left out variable, otherwise need to create a new factor
  # if(length(anyleftout) == 1){
  #   crossloadmodel_clean[[num_factor]] <- paste0(crossloadmodel_clean[[num_factor]], '+', anyleftout)
  # }
  if(length(anyleftout) == 1){
    crossloadmodel_clean[[num_factor]] <- paste0(crossloadmodel_clean[[num_factor]], '+', anyleftout)
  }
  if(length(anyleftout) == 0){
    crossloadmodel_clean <- crossloadmodel_clean
  }
  if(length(anyleftout)>1){
    #setup stepPrev object to select scaling indiactor for this new factor
    stepPrev <- list()
    stepPrev$goodmodelpart <- crossloadmodel_clean
    stepPrev$badvar <- anyleftout
    stepPrev$num_factor <- num_factor
    stepPrev$correlatedErrors <- correlatedErrors
    scalingind <- select_scalingind_stepN(data, sigLevel, scalingCrit,
                                          stepPrev = stepPrev)
    crossloadmodel_clean[[num_factor+1]] <- paste0('f', num_factor+1, '=~',
                                                   scalingind, '+',
                                                   paste0(anyleftout[-which(anyleftout==scalingind)], collapse = '+'))
  }


  ##new fit
  newcrossloadfit <- miive(model = paste0(paste0(crossloadmodel_clean, collapse = '\n'), '\n', correlatedErrors), data, var.cov = T)

  ##get badvar
  newbadvar <- getbadvar(newcrossloadfit, sigLevel = sigLevel)


  ##decide if need next step
  nextstep <- 'no'
  # if(length(newbadvar) > 1 | length(anyleftout) > 1){
  #   nextstep <- 'yes'
  # }
  if(length(newbadvar) > 1 && !identical(newbadvar,badvar)){
    nextstep <- 'yes'
    finalmodel <- paste0(crossloadmodel_clean, collapse = '\n')
    ##create the final obj to return
    finalobj <- list(model = finalmodel,
                     fit = newcrossloadfit,
                     num_factor = num_factor,
                     badvar = c(newbadvar, anyleftout), #need to combine the badvars and left out variables
                     nextstep = nextstep,
                     correlatedErrors = correlatedErrors,
                     model_sep = crossloadmodel_sepclean,
                     goodmodelpart = crossloadmodel_clean)
  }
  if(identical(newbadvar,badvar)){ #if new model still has the same bad variables as the previous step
    #then we go back to the previous step and stop the algorithm
    finalobj <- stepPrev
    finalobj$nextstep <- 'no'
  # }else{
  #   ##include correlated errors when applicable for final model print
  #   if(!is.null(correlatedErrors)){
  #     finalmodel <- paste0(paste0(crossloadmodel_clean, collapse = '\n'),
  #                          '\n', correlatedErrors)
  #   }else{finalmodel <- paste0(crossloadmodel_clean, collapse = '\n')}
    # finalmodel <- paste0(crossloadmodel_clean, collapse = '\n')
    # ##create the final obj to return
    # finalobj <- list(model = finalmodel,
    #                  fit = newcrossloadfit,
    #                  num_factor = num_factor,
    #                  badvar = c(newbadvar, anyleftout), #need to combine the badvars and left out variables
    #                  nextstep = nextstep,
    #                  correlatedErrors = correlatedErrors,
    #                  model_sep = crossloadmodel_sepclean,
    #                  goodmodelpart = crossloadmodel_clean)
  }


  return(finalobj)

}

EFAmiiv_multi <- function(model = NULL,
                          data,
                          sigLevel = .05,
                          scalingCrit = "sargan+factorloading_R2",
                          correlatedErrors = NULL){

  s1 <- EFAmiiv_multi_s1(model = model,
                          data = data,
                          sigLevel = sigLevel,
                          scalingCrit = scalingCrit,
                          correlatedErrors = correlatedErrors)

  if(s1$nextstep == 'no'){
    finalobj <- s1
  }
  if(s1$nextstep == 'yes'){
    s2 <- EFAmiiv_multi_s2(s1, data, sigLevel, scalingCrit)
    finalobj <- s2
    while (s2$nextstep == 'yes') {
      s2 <- EFAmiiv_multi_s2(s1, data, sigLevel, scalingCrit)
      finalobj <- s2
    }
  }
  if(!is.null(finalobj$correlatedErrors)){
    finalmodel <- paste0(model, '\n', paste0(finalobj$correlatedErrors, collapse = '\n'))
  }
  if(is.null(finalobj$correlatedErrors)){
    finalmodel <- finalobj$model
  }

  return(list(model = finalmodel,
              fit = finalobj$fit))

}
