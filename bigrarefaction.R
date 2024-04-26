# Calculations for interpolation with big numbers

# (c) Oleksii Dubovyk, oadubovyk@gmail.com
# Dept of Biological Sciences, Old Dominion University, Norfolk, VA 23529, USA
# Apr 26, 2024

library(gmp)

interpolation <- function(abundances, ind, N = sum(abundances), S = length(abundances), mode = "basic"){
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Analytical estimation of interpolated species richness
  #
  # Args:
  # abundances - vector of nonzero species abundances
  # ind - target number of individuals
  # N - optional, total number of individuals in the community
  # S - optional, total number of species in the community
  # mode - calculation mode:
  #   "base" or "b" - classic calculation with binomial coefficients
  #   "alternative" or "a" - alternative formula avoiding binomial coefficients but still generating large numbers
  #   "large" or "l" - alternative approach invoking gmp package
  #
  # Output:
  # Numeric, expected species richness.
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if (ind < 0) stop("interpolation :: target number of individuals is < 0")
  if (ind > N) stop("interpolation :: target number of individuals > observed, use extrapolation instead")
  S_exp <- 0
  if (mode == "basic" | mode == "b"){
    S_exp <- S - choose(N, ind)^(-1) * sum(sapply(abundances, function(i) choose(N - i, ind)))
  }else if (mode == "alternative" | mode == "a"){
    S_exp <- S - sum(sapply(abundances, function(i) prod((N-i-ind+1):(N-ind))/prod((N-i+1):(N))))
  }else if (mode == "large" | mode == "l"){
    S_exp <- sapply(abundances, function(i) gmp::prod.bigz((N-i-ind+1):(N-ind))/gmp::prod.bigz((N-i+1):(N))) %>% gmp::c_bigq() %>% gmp::sum.bigq()
    S_exp <- S - as.numeric(S_exp)
  }
  return(S_exp)
}

extrapolation <- function(abundances, 
                          ind, 
                          N = sum(abundances), 
                          S = length(abundances), 
                          f0 = (length(abundances[abundances == 1])^2) / (2*length(abundances[abundances == 2])), 
                          f1 = length(abundances[abundances == 1])){
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Analytical estimation of extrapolated species richness through Chao estimator
  #
  # Args:
  # abundances - vector of nonzero species abundances
  # ind - target number of individuals
  # N - optional, total number of individuals in the community
  # S - optional, total number of species in the community
  # f0 - optional, Chao1 estimator
  # f1 - optional, number of singleton species
  #
  # Output:
  # Numeric, expected species richness.
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if (ind < 0) stop("extrapolation :: target number of individuals is < 0")
  if (ind < N) stop("extrapolation :: target number of individuals < observed, use interpolation instead")
  m <- ind - N
  S_exp <- S + f0*(1-(1-(f1/(N*f0 + f1)))^m)
  return(S_exp)
}

rarefaction <- function(abundances, 
                        ind, 
                        N = sum(abundances), 
                        S = length(abundances), 
                        f0 = (length(abundances[abundances == 1])^2) / (2*length(abundances[abundances == 2])), 
                        f1 = length(abundances[abundances == 1])){
  
  if (ind <= 0){
    return(0)
  }else{
    S_exp <- numeric(1)
    if (ind < N){
      S_exp <- interpolation(abundances = abundances, ind = ind, N = N, S = S, mode = "l")
    }else if (ind == N){
      S_exp = S
    }else{
      m <- ind - N
      S_exp <- extrapolation(abundances = abundances, ind = ind, N = N, S = S, f0 = f0, f1 = f1)
    }
    return(S_exp)
  }
}