# === The objective here is to implement the k_samples testing procedure in Chen and Hansen

# There are two algorithms to specify
# Firstly, the one to obtain the marginal density p(x \ c ) where c is the total mass of the central measure
# p(x \ c) = normalpdf * q_theta(x\c)

# Implement algorithm 1 in the paper where
# J = depth of the Polya tree
# c =  total mass of the prior central measure
# x = vector of data
# also have mu and sigma as parameters of the normal central distribution
# The algorithm returns log(q) 

algo1 <- function(J, c, x, log = TRUE, mu, sigma){
  
  # mu = mean(x)
  # sigma = sqrt( sum((x - mu)^2 )/length(x))
  
  s = lgamma(2*c) - lgamma(2*c + length(x)) - 2^J * lgamma(c*J^2)
  
  # count in which partition observation fall
  n_k <- rep(0, 2^J)
  
  for (i in 1:length(x)){
    n_k[ ceiling(2^J*pnorm(x[i], mean = mu, sd = sigma))] = n_k[ ceiling(2^J*pnorm(x[i], mean = mu, sd = sigma))] + 1
  }
  
  for (nk in n_k){
    s = s + lgamma(c*J^2 + nk)
  }
  
  if (J > 1){
    for (j in (J-1):1){
      n_k <- n_k + head(c(0, n_k),-1)
      n_k <- n_k[seq_along(n_k)%%2 == 0] 
      
      for (nk in n_k){
        # this gives NaN too easily because it doesn t simplify
        # s = s + lgamma(2*c*(j + 1)^2) + lgamma(c*j^2 + nk) - lgamma(c*j^2) - lgamma(2*c*(j+1)^2+nk)
        
        if(nk != 0){
          s = s + sum(log(c*j^2 + 0:(nk-1))) - sum(log(2*c*(j+1)^2 + 0:(nk-1)))
        }
        
      }
      
      }
  }
  
  if(log){return(s)
    }else{return(exp(s))}
}


# == We can now focus on computing Bayes Factor.
# We need a function that computes the log numerator. This will need to be used over and over for different permutations 
# And need a function that computes the log denominator (which only needs to be computed once)

log_Tobs_numerator <- function(x1, x2, J){
  
  # compute MLE estimators for the centering dbution for each sample
  mu1 = mean(x1); mu2 = mean(x2)
  sd1 = sqrt( sum((x1 - mu1)^2)/length(x1))
  sd2 = sqrt( sum((x2 - mu2)^2)/length(x2))
  
  # Set up maximisation of c1 and c2 (mass of centering dbution)
  c1 = c2 = 1
  M1 = M2 = -Inf
  
  # maximise q(x|c) over grid: exp 19 (j − 1) − 7 
  for (i in (1:20)){
    tmp1 <- algo1(J = J, c = exp((14*((i-1)-7)/19)), x = x1)
    tmp2 <- algo1(J = J, c = exp((14*((i-1)-7)/19)), x = x2)

    if (tmp1 > M1){ c1 = exp((14*((i-1)-7)/19)); M1 = tmp1}
    if (tmp2 > M2){ c2 = exp((14*((i-1)-7)/19)); M2 = tmp2}
  }
  
  # log of each component 
  out1 <- length(x1)*J*log(2) + sum(dnorm(x1, mean = mu1, sd = sd1, log = TRUE)) + M1
  out2 <- length(x2)*J*log(2) + sum(dnorm(x2, mean = mu2, sd = sd2, log = TRUE)) + M2
  
  return(out1 + out2)
}


log_Tobs_denominator <- function(x0, J){
  
  # basically same as above but we only need to do it once (permutation invariance)
  mu0 = mean(x0)
  sd0 = sqrt( (x0-mu0)^2/length(x0) )
  
  c0 = 1; M0 = -Inf
  
  for (i in (1:20)){
    tmp0 <- algo1(J = J, c = exp((14*((i-1)-7)/19)), x = x0)
    if (tmp0 > M0){ c0 = exp((14*((i-1)-7)/19)); M0 = tmp0}
  }
  
  return( length(x0)*J*log(2) + sum(dnorm(x0, mean = mu0, sd = sd0, log = TRUE)) + M0 )
  
}







permutation_test <- function(x1, x2, J, M = 1000 ){
  
  x0 <- c(x1, x2)
  N0 = length(x0); N1 = length(x1)
  
  # Get the Test statistic (Do we even need to calculate the denominator)
  log_Tobs <- log_Tobs_numerator(x1, x2, J) - log_Tobs_denominator(x0, J)
  
  # Approximate distribution of statistic under null hypothesis
  test_dbution <- rep(NA, M)
  for (i in 1:M){
    
    indices <- sample(1:N0)[1:N1]
    x1_tmp <- x0[indices]
    x2_tmp <- x0[-indices]
    test_dbution[i] <- log_Tobs_numerator(x1_tmp, x2_tmp, J)
  }
  test_dbution <- test_dbution - log_Tobs_denominator(x0, J)
  return(mean(log_Tobs > test_dbution))
  
}




# https://rdrr.io/cran/NADA/man/cenmle.html


