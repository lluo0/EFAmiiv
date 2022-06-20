
getparam <- function(model,
                     var_obs = 1.25, #observed variables' error variance set to be 1.25 by default
                     var_error = NULL,
                     R2 = NULL){
  latentFac <- observedVar <- vector()
  ######setup#####
  modelstring <- strsplit(model, "\n")[[1]]
  #get the number of equations specifying loadings. no correlated errors here.
  latentEq <- as.vector(which(sapply(modelstring, function(u) grepl('=~', u))))
  #get the names of the latent factors. !use gsub to get rid of spaces
  latentFac <- gsub(' ','',as.vector(sapply(modelstring[latentEq], function(i) strsplit(i, '=~')[[1]][1])))
  # #the number of equations with correlated errors for latent factors
  # corrEq_latentFac <- as.vector(which(sapply(modelstring,
  #                                            function(u) grepl(paste(latentFac, collapse = '|'), u))))[-latentEq]
  # #the number of equations with correlated errors for observed variables
  # corrEq_observedVar <- c(1:length(modelstring))[-c(latentEq, corrEq_latentFac)]

  #####get names and loadings for observed variables#####
  observedVar_list <- lapply(modelstring[latentEq],
                             function(u) gsub(".*[=~]([^.]+)", "\\1", u)) #get strings after =~
  observedVar_list <- lapply(observedVar_list,
                             function(u) strsplit(u, '+', fixed = T)[[1]]) #separate strings by +

  observedVar <- lapply(observedVar_list,
                        function(u) as.vector(sapply(u,
                                                     function(i) gsub(' ', '', sapply(i, #get rid of spaces
                                                                                      function(m) gsub(".*[*]([^.]+)", "\\1", m)))))) #get strings after *
  observedVar_loadings <- lapply(observedVar_list,
                                 function(u) as.numeric(sapply(u,
                                                               function(i) gsub(' ', '', sapply(i, #get rid of spaces
                                                                                                function(m) gsub("[*].*$", "", m)))))) #get strings before *
  observedVar_unique <- unique(unlist(observedVar))
  ######create the loading matrix#####
  Lambda <- matrix(0,
                   nrow = length(observedVar_unique),
                   ncol = length(latentFac))
  rownames(Lambda) <- observedVar_unique #note that the order of observed variables is the same as they appear in the model
  colnames(Lambda) <- latentFac
  for(n in 1:ncol(Lambda)){
    for(p in 1:nrow(Lambda)){
      if(rownames(Lambda)[p] %in% observedVar[[n]]){
        ordernum <- which(observedVar[[n]] == rownames(Lambda)[p]) #get the corresponding loadings' order in that list
        Lambda[p,n] <- as.numeric(observedVar_loadings[[n]][ordernum])
      }
    }
  }

  #####get the equation numbers for var-cov matrices#####
  #the number of equations with correlated errors for observed variables
  corrEq_observedVar <- as.vector(which(sapply(modelstring,
                                               function(u) grepl(paste(observedVar_unique, collapse = '|'), u))))[-latentEq]
  #the number of equations with correlated errors for latent factors
  corrEq_latentFac <- c(1:length(modelstring))[-c(latentEq, corrEq_observedVar)]

  ######get covariance of latent factors######
  varcov_latent <- diag(1, nrow = length(latentFac)) #variance of latent factors set to be 1 and 0 covariance if not specified otherwise
  colnames(varcov_latent) <- rownames(varcov_latent) <- latentFac
  # if(length(corrEq_latentFac)==0){ #no correlation between latent factors
  #   varcov_latent <- matrix(0, nrow = length(latentFac), ncol = length(latentFac))
  #   colnames(varcov_latent) <- rownames(varcov_latent) <- latentFac
  # }
  if(length(corrEq_latentFac)!=0){ #get the correlation between latent factors based on specified model
    latent_list <- lapply(modelstring[corrEq_latentFac],
                          function(u) strsplit(u, '~~', fixed = T)[[1]]) #separate strings by ~~
    for(n in 1:length(latent_list)){
      names(latent_list)[n] <- gsub(".*[*]([^.]+)", "\\1", latent_list[[n]][2]) #name each list the other latent factor
      latent_list[[n]][1] <- gsub(' ', '', latent_list[[n]][1]) #get rid of spaces
      latent_list[[n]][2] <- as.numeric(gsub("[*].*$", "", latent_list[[n]][2])) #leave the 2nd element the covariance
    }
    ######fill the latent factor var-cov matrix#####
    for(n in 1:length(latent_list)){
      for(p in 1:nrow(varcov_latent)){
        if(names(latent_list)[n]==rownames(varcov_latent)[p]){
          for(i in 1:ncol(varcov_latent)){
            # varcov_latent[p,i] <- varcov_latent[i,p] <- ifelse(colnames(varcov_latent)[i]==latent_list[[n]][1],
            #                              as.numeric(latent_list[[n]][2]), 0)
            varcov_latent[i,p] <-  ifelse(colnames(varcov_latent)[i]==latent_list[[n]][1],
                                          as.numeric(latent_list[[n]][2]), varcov_latent[i,p])
            varcov_latent[p,i] <-  ifelse(colnames(varcov_latent)[i]==latent_list[[n]][1],
                                          as.numeric(latent_list[[n]][2]), varcov_latent[p,i])
          }
        }
      }
    }
  }
  ######get covariance of observed variables errors#####
  varcov_epsilon <- matrix(0, nrow = length(observedVar_unique), ncol = length(observedVar_unique) )
  colnames(varcov_epsilon) <- rownames(varcov_epsilon) <- observedVar_unique
  # #print error when none of error covariance, observed variable variance, or R2 is specified
  # if(is.null(var_error) && is.null(var_obs) && is.null(R2)) stop('Please specify at least one of the following: \n
  #                                                                var_error: observed variables error variance, \n
  #                                                                var_obs: observed variables variance, or \n
  #                                                                R2: r squared values for observed variables.')

  if(is.null(var_error) && is.null(R2)){ #calculate error variance when no other information specified using default variance of observed variables of 1.25
    var_latent <- diag(varcov_latent)
    var_error <- var_obs - as.vector(Lambda^2%*%as.matrix(var_latent)) #V(e) = V(x) - lambda^2*V(latent)
    diag(varcov_epsilon) <- var_error
  }
  if(!is.null(var_error)){ #include variance of error terms when specified directly
    diag(varcov_epsilon) <- var_error
  }
  if(!is.null(R2)){ #calculate error variance when R2 specified
    var_latent <- diag(varcov_latent)
    var_obs <- as.vector(Lambda^2%*%as.matrix(var_latent))/R2 #V(x) = lambda^2*V(latent)/R2
    var_error <- var_obs - as.vector(Lambda^2%*%as.matrix(var_latent)) #V(e) = V(x) - lambda^2*V(latent)
    ##QUESTION: does the equation changes when there's correlated errors present???
    diag(varcov_epsilon) <- var_error
  }

  # if(length(corrEq_observedVar)==0){
  #   varcov_epsilon <- matrix(0, nrow = length(observedVar_unique), ncol = length(observedVar_unique) )
  #   colnames(varcov_epsilon) <- rownames(varcov_epsilon) <- observedVar_unique
  # }
  if(length(corrEq_observedVar)!=0){ #when correlated errors for observed variables present, need to get off-diagonal values
    epsilon_list <- lapply(modelstring[corrEq_observedVar],
                           function(u) strsplit(u, '~~', fixed = T)[[1]]) #separate strings by +
    for(n in 1:length(epsilon_list)){
      names(epsilon_list)[n] <- gsub(' ','',gsub(".*[*]([^.]+)", "\\1", epsilon_list[[n]][2])) #name each list the other latent factor
      epsilon_list[[n]][1] <- gsub(' ', '', epsilon_list[[n]][1]) #get rid of spaces
      epsilon_list[[n]][2] <- as.numeric(gsub("[*].*$", "", epsilon_list[[n]][2])) #leave the 2nd element the covariance
    }
    # ######create latent factor var-cov matrix#####
    # varcov_epsilon <- diag(0, nrow = length(observedVar_unique))
    # colnames(varcov_epsilon) <- rownames(varcov_epsilon) <- observedVar_unique
    #fill in var-cov values from the model specification
    for(n in 1:length(epsilon_list)){
      for(p in 1:nrow(varcov_epsilon)){
        if(names(epsilon_list)[n]==rownames(varcov_epsilon)[p]){
          for(i in 1:ncol(varcov_epsilon)){
            varcov_epsilon[p,i] <- ifelse(colnames(varcov_epsilon)[i]==epsilon_list[[n]][1],
                                          as.numeric(epsilon_list[[n]][2]), varcov_epsilon[p,i])
            varcov_epsilon[i,p] <- ifelse(colnames(varcov_epsilon)[i]==epsilon_list[[n]][1],
                                          as.numeric(epsilon_list[[n]][2]), varcov_epsilon[i,p])
          }
        }
      }
    }
  }
  #
  # # diag(varcov_epsilon) <- 1-rowSums(Lambda^2) #standardize
  # diag(varcov_epsilon) <- epsilon
  finalobj <- list(Lambda = Lambda,
                   varcov_latent = varcov_latent,
                   varcov_epsilon = varcov_epsilon)
  return(finalobj)
}

########examples######
model <- 'f1=~1*x1 + .8*x2 + .7*x3 + .7*x4 + .3*x7
f2=~ 1*x5 + .8*x6 + .7*x7 + .7*x8
f3=~ 1*x9 + .8*x10 + .7*x11 + .6*x12
f1 ~~ .5*f2
f1 ~~ .4*f3
f2 ~~ .6*f3
x2 ~~ .3*x3'
getparam(model, R2 = rep(.8, 12))
