# Want to run specific experiment and accept some arguments

require(optparse)

option_list <- list(
  make_option("--rate", default = 0.5,
              help = 'Experimental Censoring Rate')
)

opt <- parse_args(OptionParser(option_list=option_list))


# load stuff to perform HMC permutation test:

source('k_samples_testing.R')
source('censored_MCMC.R')
source('censored_permutation_test.R')


# Here we will need to load/simulate data


# Then we will need to apply our HMC test as well as competitors




# == What type of experiments do we want???
# Here I compare non censored with censored


f <- function(x){
  1-exp(-x)-x*opt$rate
}
M <- rootSolve::uniroot.all(f, interval = c(0.01,100))

results <- data.frame(iteration = as.numeric(),
                      p_value = as.numeric())

for (i in 1:100){
  
  print(paste0('Simulation number: ', i))
  
  # simulate uncensored data 
  t1 <- rexp(50)
  delta1 <- rep(1, 50)
  
  # simulate censored data
  x2 <- rexp(50)
  c2 <- runif(50, min = 0, max = M)
  t2 <- mapply(min,x2, c2)
  delta2 <- as.numeric( x2 < c2 )
  
  # perform permutation test
  p <- censored_permutation_test(t1 = t1, t2 = t2, delta1 = delta1, delta2 = delta2)$p_value
  results <- rbind(results, c(i, p))
  
  
}

results$censoring_rate <- opt$rate

opt$rate <- gsub('\\.', '_', opt$rate)

filename <- paste0('type1error_MH_rate_', opt$rate, '.rds')
saveRDS(results,file=filename)


