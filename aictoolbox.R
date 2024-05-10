metaAIC <- function(model){
  #
  # get AIC parameters
  #
  aic <- model$aic %>% unlist()
  lL <- logLik(model) %>% unlist()
  df <- sum(model$edf) %>% unlist()
  f <- model$formula %>% unlist()
  return(c("formula" = f, "AIC" = aic, "logLik" = lL, "df" = df))
}

rankAIC <- function(models){
  aicweights <- function(aics){
    daic <- aics - min(aics)
    w <- numeric(length(aics))
    s <- sum(exp(-0.5*daic))
    for (i in 1:length(aics)){
      w[i] <- exp(-0.5*daic[i])/s
    }
    return(w)
  }
  
  formulae <- character(0)
  aics <- numeric(0)
  lLs <- numeric(0)
  dfs <- numeric(0)
  for (model in models){
    formulae <- c(formulae, model$formula %>% paste(collapse = " "))
    aics <- c(aics, model$aic)
    lLs <- c(lLs, scam::logLik.scam(model))
    dfs <- c(dfs, sum(model$edf))
  }
  aictable <- tibble(model = formulae,
                     AIC = aics,
                     logLik = lLs,
                     df = dfs)
  aictable$delta_AIC <- aictable$AIC - min(aictable$AIC)
  aictable$weight <- aicweights(aictable$AIC)
  aictable[,c("model", "AIC", "delta_AIC", "weight", "logLik", "df")] %>%
    as_tibble() %>%
    arrange(delta_AIC)
}
