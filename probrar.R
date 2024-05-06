# Probabilistic rarefaction

# (c) Oleksii Dubovyk, oadubovyk@gmail.com
# Dept of Biological Sciences, Old Dominion University, Norfolk, VA 23529, USA
# May 05, 2024

library(tidyverse)
library(gmp)

prob_same <- function(N, m){
  S <- length(N)
  J <- sum(N)
  P <- sapply(N, function(x) x/J)
  N <- P*J
  sum(sapply(1:S, function(i) prod(sapply(1:m, function(k){
    (N[i]-P[i]*(k-1))/(J-(k-1))
  }))))
}

probs_roll <- function(N){
  sapply(1:sum(N), function(m){
    prob_same(N = N, m = m)
  })
}

probs_combin <- function(n, k = 1){
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Probability of encountering a new species by drawing a new individual
  #
  # Args:
  # n - numeric vector, abundances of species in a community
  # k - step size in drawing individuals
  #
  # Output:
  # Numeric, estimated probabilities of encountering a new species at [i]th individual.
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  test1 <- function(n, m){
    N <- sum(n)
    S <- length(n)
    sapply(1:S, function(i) -((gmp::chooseZ(N-n[i], m))/(gmp::chooseZ(N, m)))+((gmp::chooseZ(N-n[i], m-1))/(gmp::chooseZ(N, m-1))))
  }
  if(k == 1){
    sapply(1:sum(n), function(x) test1(n, x) %>% gmp::c_bigq() %>% gmp::sum.bigq() %>% as.numeric())
  }else{
    l <-  floor(sum(n)/k)
    idx <-  as.data.frame(sapply(1:l, function(x) seq(x*k-(k-1), x*k, 1)))
    probs <- as.data.frame(apply(idx, 2, function(x) sapply(x, function(y) sum(test1(n=n, m=y)))))
    lvls <- nrow(probs)
    cum_probs <- numeric(l)
    for (p in 1:ncol(probs)){
      formula <- c()
      for (i in 1:lvls){
        formula <- c(formula, 
                     ifelse(i%%2==1, 1, -1)*sum(sapply(1:nrow(expand.grid(lvls, i)), function(j) prod(sapply(expand.grid(lvls, i)[j,], function(k) probs[k, p])))))
      }
      cum_probs[p] <- sum(formula)
    }
    cum_probs
  }
}

beyond_combin <- function(P, ceil, burnin = round(0.5*length(P)), threshold = 1e-30){
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Probabilistic extrapolation
  #
  # Args:
  # P - estimated probabilities of encountering a new species at [i]th individual,
  #   an output of probs_combin()
  # ceil - ceiling until which to continue extrapolation (as a number of drawn individuals).
  #   If missing, the extrapolation will continue until there is no further increase
  #   in species richness within threshold specified by `threshold` argument.
  # burnin - burn-in steps until which the first elements of `P` should be ignored when
  #   estimating how probabilities decrease with the sample size.
  # threshold - threshold by which it is assumed that there is no further increase in 
  #   probabilities estimates. Setting to 0 is not advised since the probabilities can
  #   exhibit very small fractions.
  #
  # Output:
  # Numeric, estimated probabilities of encountering a new species at [i]th individual.
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  probs_obs <- P
  
  if (missing(ceil)){
    
    if (probs_obs[length(probs_obs)] <= threshold){
      out <- probs_obs[probs_obs > 0]
    }else{
      probs_considered <- probs_obs[burnin:length(probs_obs)]
      factor1 <- sapply(2:length(probs_considered), function(i) probs_considered[i]/probs_considered[i-1])
      factor2 <- sapply(2:length(factor1), function(i) factor1[i]/factor1[i-1])
      factor3 <- sapply(2:length(factor2), function(i) factor2[i]/factor2[i-1])
      out <- probs_obs
      out[length(out) + 1] <- out[length(out)] * mean(factor1) * mean(factor2) * mean(factor3)
      p_i <- probs_obs[length(probs_obs)]
      while (p_i > threshold){
        p_i <- out[length(out)] * mean(factor1) * mean(factor2) * mean(factor3)
        out <- c(out, p_i)
      }
    }
    
  }else{
   
    if (ceil <= 0){
      
      out <- 0
      
    }else if (ceil <= length(P)){
      
      out <- probs_obs[0:ceil]
      
    }else{
      
      probs_obs <- probs_obs[probs_obs > 0]
      probs_considered <- probs_obs[burnin:length(probs_obs)]
      factor1 <- sapply(2:length(probs_considered), function(i) probs_considered[i]/probs_considered[i-1])
      factor2 <- sapply(2:length(factor1), function(i) factor1[i]/factor1[i-1])
      factor3 <- sapply(2:length(factor2), function(i) factor2[i]/factor2[i-1])
      # very brute force to approximate a series
      out <- probs_obs
      out[length(out) + 1] <- out[length(out)] * mean(factor1) * mean(factor2) * mean(factor3)
      for (i in (length(out)+1):ceil){
        out[i] <- out[i-1] * mean(factor1) * mean(factor2) * mean(factor3)
      }
      
    }
     
  }
  
  return(out)
  
}

beyond_combin(probs_combin(c(4, 4))) %>% cumsum() %>% plot(type = "l")

# # Example:
# test_obs <- probs_combin(c(4, 3, 2, 1))
# beyond_combin(test_obs)
