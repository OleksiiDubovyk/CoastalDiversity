# Probabilistic rarefaction

# (c) Oleksii Dubovyk, oadubovyk@gmail.com
# Dept of Biological Sciences, Old Dominion University, Norfolk, VA 23529, USA
# Apr 26, 2024

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

beyond_combin <- function(N, ceil = round(2*sum(N)), burnin = round(0.5*sum(N))){
  probs_obs <- probs_combin(N)
  
  if (ceil <= 0){
    0
  }else if (ceil <= sum(N)){
    probs_obs[0:ceil]
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
    return(out)
  }
}
