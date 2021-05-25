# Implement the censored permutation test 

require(parmsurvfit)
source('censored_MCMC.R')
source('k_samples_testing.R')

check_arguments <- function(t, delta){
  
  if (length(t) != length(delta)){
    stop("Lengths of delta and t should match")
  }
  
  if (sum(!(delta == 1 | delta ==0)) != 0){
    stop("delta should be a binary vector")
  }
}


censored_permutation_test <- function(t1, t2, delta1, delta2, N_permutations = 200){
  
  # simple check
  check_arguments(t1, delta1)
  check_arguments(t2, delta2)
  
  # null hypothsis dataset
  t0 <- c(t1, t2)
  delta0 <- c(delta1, delta2)
  
  L1 <- length(t1)
  L0 <- length(t0)
  J <- ceiling( log(L0, base = 2) )
  
  
  # compute \bar{p_i} for censored groups and datasets
  # bar_p0 <- censored_MH(t = t0, delta = delta0, J = J)$Ep # Do we actually need this??
  bar_p1 <- censored_MH(t = t1, delta = delta1, J = J)$Ep
  bar_p2 <- censored_MH(t = t2, delta = delta2, J = J)$Ep
  
  statistic_numerator <- bar_p1 * bar_p2
  
  
  numerator_distribution <- rep(0, N_permutations)
  # now shuffle
  for (i_perm in 1:N_permutations){
    
    # print(paste0('Permutation Number: ', i_perm, '\r'))
    
    # permute and divide in new groups with same sizes as t1, t2
    new_groups_indices <- sample(1:L0)
    group1_indices <- new_groups_indices[1:L1]
    group2_indices <- new_groups_indices[(L1+1):L0]
    
    new_t1 <- t0[group1_indices]
    new_t2 <- t0[group2_indices]
    
    new_delta1 <-delta0[group1_indices]
    new_delta2 <-delta0[group2_indices]
    
    # compute numerator of stat for this partition
    new_bar_p1 <- censored_MH(t = new_t1, delta = new_delta1, J = J)$Ep
    new_bar_p2 <- censored_MH(t = new_t2, delta = new_delta2, J = J)$Ep
    
    # store statistic numerator
    numerator_distribution[i_perm] <- new_bar_p1 * new_bar_p2
  }
  
  p_value <- mean(statistic_numerator < numerator_distribution )
  
  return(list(p_value = p_value, statistic_numerator = statistic_numerator, numerator_distribution = numerator_distribution))
}


