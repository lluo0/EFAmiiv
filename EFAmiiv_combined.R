

EFAmiiv_cbnd <- function(model = NULL,
                         data,
                         sigLevel = .05,
                         scalingCrit = "sargan+factorloading_R2",
                         correlatedErrors = NULL){
  if(is.null(model)){
    finalobj <- EFAmiiv(data = data,
                        sigLevel = sigLevel,
                        scalingCrit = scalingCrit,
                        correlatedErrors = correlatedErrors)
  }
  if(!is.null(model)){
    finalobj <- EFAmiiv_multi(model = model,
                              data = data,
                              sigLevel = sigLevel,
                              scalingCrit = scalingCrit,
                              correlatedErrors = correlatedErrors)
  }
  return(finalobj)
}
