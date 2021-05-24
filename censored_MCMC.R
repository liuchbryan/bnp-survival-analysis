# == Objective of this script is to implement algorithm 2 in Chen and Hanson


# == Data ==  

# There are I different groups of observations, each with n_i datapoints
# Observations for individual j in group i is (^x_ij, delta_ij).
# x will denote the vector of logged deaths' times (some obs will be latent)

# If observation is uncensored (delta_ij = 1):
# x_ij <- log(t_ij)
# Else:
# x_ij \in (log(t_ij), \infty)

# Need to be able to perform MLE estimates of Normal distribution censored samples,
# look at:  https://rdrr.io/cran/NADA/man/cenmle.html


# unfeasible to compute the statistic as before (as for any observation of latent vbles would need permutation test)
# Instead consider:
# E [ p_theta (x_i \ c_i) \ D_i]
# Also, computationally unfeasible to maximise this expectation wrt c, so put a prior instead.
# c_i ~ Gamma(5,1)





#== HELPER FUNCTIONS == # --------------------


require('parmsurvfit')
source('k_samples_testing.R')

sample_censored <- function(x,censored_indices, mu_hat, sd_hat){
  
  a <- x; b <- x
  b[censored_indices] <- Inf 
  
  censored_length <- length(censored_indices)
  
  u <- runif(censored_length, min = pnorm(a, mean = mu_hat, sd = sd_hat),
             max = pnorm(b, mean = mu_hat, sd = sd_hat)) # could just set this to one
  
  x[censored_indices] <- qnorm(u, mean = mu_hat, sd = sd_hat)
  
  return(x)
}

# function that takes x, x_prop -> determines s0
compute_s0 <- function(xj, xj_prop, J, mu_hat, sd_hat){
  
  for (level in 1:J){
    
    bin1 <- ceiling(2^level*pnorm(xj, mean = mu_hat, sd = sd_hat))
    bin2 <- ceiling(2^level*pnorm(xj_prop, mean = mu_hat, sd = sd_hat))
    
    if (bin1 != bin2){return(level-1)}
  }
  return(J)
}

# --------------------------------------------

censored_MH <- function(t, delta, J, B = 100, M = 200, V = 1.5){
  
  
  # check lengths
  N <- length(t)
  if (length(delta) != N){stop("Lengths of delta and t should match")}
  
  # Transform to log scale and define bounds for each observation
  censored_indices <- which(delta == 0)
  censored_length <- length(censored_indices) 
    
  x <- log(t)
  
  # compute MLEs 
  # do i need to do it here?
  
  data <- data.frame(time = t, censor = delta)
  fit <- fit_data(data, dist = 'lnorm')
  
  mu_hat <- fit$estimate[1]
  sd_hat <- fit$estimate[2]
  
  # sample latent log censored deaths
  x <- sample_censored(x, censored_indices, mu_hat, sd_hat)
  
  # initialise stuff
  c <- 1
  SUM <- 0
  
  # Compute n(J, k) matrix of 'bin' counts of the 'normal' partition
  n <- list() 
  n[[J]] <- rep(0, 2^J)
  for (i in 1:N){
    bin <- ceiling(2^J*pnorm(x[i], mean = mu_hat, sd = sd_hat))
    n[[J]][bin] <- n[[J]][bin] + 1
  }
  
  if (J > 1){
    for (j in (J-1):1){
      n[[j]] <- rep(0, 2^j)
      
      n[[j]] <- n[[j+1]] + head(c(0, n[[j+1]]),-1)
      n[[j]] <- n[[j]][seq_along(n[[j]])%%2 == 0] 
    }
  }
  
  # Helper function
  
  # MH
  
  log_qx <- algo1(J = J, c = c, x = x, mu = mu_hat, sigma = sd_hat)
  
  for( iter in 1:(B+M)){
    
    # STEP 1 # -------------------------------------------
    
    x_prop <- sample_censored(x, censored_indices, mu_hat, sd_hat)
    # pwise comparison x_prob[censored_indices] vs x[censored_indices]
    
    k <- ceiling(2^J*pnorm(x, mean = mu_hat, sd = sd_hat))[censored_indices]
    k_prop <- ceiling(2^J*pnorm(x_prop, mean = mu_hat, sd = sd_hat))[censored_indices]


    s0 <- mapply(compute_s0, k, k_prop, 
           MoreArgs = list(J = J,
                           mu_hat = mu_hat,
                           sd_hat = sd_hat))
    
    
    alpha <- rep(0, censored_length)
    
    
    alpha[which(s0 != J)] <- log(c*J^2 + n[[J]][k_prop]) - log(c*J^2 + n[[J]][k])


    # What if s0 = J-1  V
    # What if s0 = J X we want the for loop not to run
    
    for (s in min(s0+1, J-1):(J-1)){
      
      # might need to make sure min(s0) < J-1 but not sure
      indices <- which(s0 < s)
      alpha[indices] <- (alpha[indices] +
                           log(c*s^2 + n[[s]][ ceiling(2^(s-J)*k_prop[indices])]) -
                           log(c*s^2 + n[[s]][ ceiling(2^(s-J)*k[indices])]) +
                           log(2*c*(s+1)^2 + n[[s]][ ceiling(2^(s-J)*k[indices])] ) -
                           log(2*c*(s+1)^2 + n[[s]][ ceiling(2^(s-J)*k_prop[indices])])
                         )
    }
    
    #(e)
    u <- runif(censored_length)
    changed_indices <- which(u < exp(alpha))
    
    # for all u 
    # x has dimension N
    # u has dimension censored_length -> ch_indices is referring to 'censored dataset' indices
    x[censored_indices[changed_indices]] <- x_prop[censored_indices[changed_indices]]

    # s0 = J 
    #bool <- ( min(s0) == J) -> might implement above
    for (s in min(s0+1, J):J){
      
      # if (bool){break}
      
      indices <- which(s0 < s)
      
      n[[s]][ceiling(2^(s-J)*k_prop[indices])] <- n[[s]][ceiling(2^(s-J)*k_prop[indices])] + 1
      n[[s]][2^(s-J)*ceiling(2^(s-J)*k[indices])] <- n[[s]][ ceiling(2^(s-J)*k[indices]) ] - 1
      
    }
  
    
    # STEP 2 etc -------------------------
    
    # compute log lik
    
    #log_px <- + sum(dnorm(x, mean = mu_hat, sd = sd_hat, log = TRUE)) + log_qx
    
    # propose c
    c_prop <- exp(rnorm(1, mean = log(c), sd = sqrt(V)))
    
    log_qx_prop <- algo1(J = J, c = c_prop, x = x, mu = mu_hat, sigma = sd_hat)
    
    log_threshold <- (log_qx_prop - log_qx +
                        dgamma(c_prop, 5,1, log = TRUE) - dgamma(c, 5,1, log = TRUE))
    
    
    
    if( log(runif(1)) < log_threshold ){
      c <- c_prop
      log_qx <- log_qx_prop
    }
    
    if (iter > B){
      log_px <- N*J*log(2) + sum(dnorm(x, mean = mu_hat, sd = sd_hat, log = TRUE)) + log_qx
      SUM <- SUM + exp(log_px)
    }
    
    
    }
  
  
  # WHAT?
  return(list(Ep = SUM/M, last_x = x))
  
  
}

# Finish it
# -> Compute MLEs  (# https://rdrr.io/cran/NADA/man/cenmle.html)
# -> Actually compute statistic and perform permutation test

2 * 200 * (300)
1 * 300



