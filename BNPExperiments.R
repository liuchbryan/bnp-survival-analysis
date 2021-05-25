library(survival)
library(IDPSurvival) # install.packages("IDPSurvival")
source('censored_permutation_test.R')

require(optparse)

option_list <- list(
  make_option("--rate", default = 0.5,
              help = 'Experimental Censoring Rate')
)

opt <- parse_args(OptionParser(option_list=option_list))

cr <- opt$rate


######################
## Experiment Stuff ##
######################

# function for simulating data from the piece-wise exponential distribution
rpwexp <- function(n, rate, intervals, cumulative = FALSE){
  if(is.null(intervals)){
    if (cumulative){return(cumsum(rexp(n,rate[1])))}else
      return(rexp(n,rate[1]))}
  k <- length(rate)
  if (k==1){
    if(cumulative){return(cumsum(rexp(n,rate)))}else
      return(rexp(n,rate))
  }
  if (length(intervals) < k-1) stop("length(intervals) must be at least length(rate) - 1")
  tx <- 0
  j <- 1
  times <- array(0,n)
  timex <- cumsum(intervals)
  indx <- array(TRUE,n)
  for(i in 1:k){
    nindx <- sum(indx)
    if (nindx==0) break
    increment <- rexp(nindx,rate[i])
    if (cumulative) times[indx] <- tx + cumsum(increment)
    else times[indx] <- tx + increment
    if (i<k){
      tx <- timex[i]
      indx <- (times > timex[i])
    }
  }
  return(times)
}

# function that implements the log-rank test 
LRTest <- function(data){
  htest <- survdiff(Surv(data$time, data$delta) ~ data$group)
  pvalue <- 1 - pchisq(q = htest$chisq, df = 1)
  decision <- ifelse(pvalue < 0.05, 1, 0)
  return(decision)
}

# function that implements the dirichlet prior test
DPTest <- function(data){
  htest <- isurvdiff(Surv(data$time, data$delta) ~ data$group, 
                     display = FALSE, nsamples = 5000, level = 0.95)
  
  decision <- htest$h
  return(decision)
}

PTTest <- function(data){
  
  # data <- data.frame(time, delta, group)
  t1 <- data$time[which(data$group == 1)]
  t2 <- data$time[which(data$group == 2)]

  delta1 <- as.integer(data$delta[which(data$group == 1)])
  delta2 <- as.integer(data$delta[which(data$group == 2)])
  
  # times, delta, group
  p_value <- censored_permutation_test(t1 = t1,t2 = t2, delta1 = delta1, delta2 = delta2)$p_value
  decision <- ifelse(p_value < 0.05, 1, 0)
  
  # Return 1 if rejects the null or 0 if accept the null
  return(decision)
}


# function that performs monte-carlo experiment for same hazards 
SHExperiment <- function(iters, n = 50, censoringrate){
  
  iter <- 0
  lrtest <- 0 # counts number of times log-rank test rejects the null
  pttest <- 0 # counts number of times polya tree prior rejects the null
  dptest <- 0 # counts number of times dirichlet process prior rejects the null
  
  lrtestind <- 0 # counts number of times log-rank test rejects the null
  pttestind <- 0 # counts number of times polya tree prior rejects the null
  dptestind <- 0 # counts number of times dirichlet process prior rejects the null
  
  while(iter < iters){
    
    print(paste0('Iteration number:', iter))
    
    group <- c(rep(1, n), rep(2, n)) # group id: treatment group id
    
    survival <- c(rexp(n, rate = 1), # survival times: treatment group 1
                  rexp(n, rate = 1)) # survival times: treatment group 2
    
    if(censoringrate == 0.00){
      time <- survival
      delta <- rep(1, n)
      data <- data.frame(time, delta, group)
      data <- data[order(data$time),]
      } else{
        if(censoringrate == 0.25){
          b = 0.606
          } else if(censoringrate == 0.50){
          b = 1.59
          } else if(censoringrate == 0.75){
          b = 3.92
          } else{
          print("Invalid Censoring Rate")
          break}
        censoring <- c(runif(n, min = 0, max = b), # censoring times: treatment group 1
                       runif(n, min = 0, max = b)) # censoring times: treatment group 2
        time <- c(NA, n + n) # time: empty vector for event times
        delta <- c(NA, n + n) # delta: survival or censoring time indicator
        
        for(i in 1:(n + n)){
          # time: event times are the first observed event 
          time[i] <- min(survival[i], censoring[i])
          
          # delta: event indicator (1: survival, 0: censoring)
          delta[i] <- ifelse(time[i] == survival[i], 1, 0)
        }
        
        # create the life table for two-sample testing
        data <- data.frame(time, delta, group)
        data <- data[order(data$time),]
      }
      
      
      # obtain decision from log-rank test and identify if it rejects the null
      lr <- LRTest(data)
      
      # obtain decision from polya-tree test and identify if it rejects the null
      pt <- PTTest(data)
      
      # obtain decision from dirichlet prior test and identify if it rejects the null
      dp <- DPTest(data)
      
      # update counters whenever we get decisive result from dirichlet process prior
      iter <- iter + 1
      
      # if dp cant decide
      if(dp == 2){
        lrtestind <- lrtestind + lr
        pttestind <- pttestind + pt
        dptestind <- dptestind + 1
      } else{
        lrtest <- lrtest + lr
        pttest <- pttest + pt
        dptest <- dptest + dp
      }
    }
    
    lrtest <- lrtest/iters
    pttest <- pttest/iters
    dptest <- dptest/iters
    
    lrtestind <- lrtestind/iters
    pttestind <- pttestind/iters
    dptestind <- dptestind/iters
    
    return(data.frame(lrtest, pttest, dptest, 
                      lrtestind, pttestind, 
                      dptestind))
}


# function that performs monte-carlo experiment for proportional hazards
PHExperiment <- function(iters, n = 50, censoringrate){
  
  iter <- 0
  lrtest <- 0 # counts number of times log-rank test rejects the null
  pttest <- 0 # counts number of times polya tree prior rejects the null
  dptest <- 0 # counts number of times dirichlet process prior rejects the null
  
  lrtestind <- 0 # counts number of times log-rank test rejects the null
  pttestind <- 0 # counts number of times polya tree prior rejects the null
  dptestind <- 0 # counts number of times dirichlet process prior rejects the null
  
  while(iter < iters){
    group <- c(rep(1, n), rep(2, n)) # group id: treatment group id
    
    survival <- c(rexp(n, rate = 1), # survival times: treatment group 1
                  rexp(n, rate = 2)) # survival times: treatment group 2
    
    if(censoringrate == 0.00){
      time <- survival
      delta <- rep(1, n)
      data <- data.frame(time, delta, group)
      data <- data[order(data$time),]
    } else{
      if(censoringrate == 0.25){
        b = c(0.606, 0.303)
      } else if(censoringrate == 0.50){
        b = c(1.59, 0.7975)
      } else if(censoringrate == 0.75){
        b = c(3.92, 1.961)
      } else{
        print("Invalid Censoring Rate")
        break
      }
      censoring <- c(runif(n, min = 0, max = b[1]), # censoring times: treatment group 1
                     runif(n, min = 0, max = b[2])) # censoring times: treatment group 2
      
      time <- c(NA, n + n) # time: empty vector for event times
      delta <- c(NA, n + n) # delta: survival or censoring time indicator
      
      for(i in 1:(n + n)){
        # time: event times are the first observed event 
        time[i] <- min(survival[i], censoring[i])
        
        # delta: event indicator (1: survival, 0: censoring)
        delta[i] <- ifelse(time[i] == survival[i], 1, 0)
      }
      
      # create the life table for two-sample testing
      data <- data.frame(time, delta, group)
      data <- data[order(data$time),]
    }
    
    # obtain decision from log-rank test and identify if it rejects the null
    lr <- LRTest(data)
    
    # obtain decision from polya-tree test and identify if it rejects the null
    pt <- PTTest(data)
    
    # obtain decision from polya-tree test and identify if it rejects the null
    dp <- DPTest(data)
    
    # update counters whenever we get decisive result from dirichlet process prior
    iter <- iter + 1
    
    if(dp == 2){
      lrtestind <- lrtestind + lr
      pttestind <- pttestind + pt
      dptestind <- dptestind + 1
    } else{
      lrtest <- lrtest + lr
      pttest <- pttest + pt
      dptest <- dptest + dp
    }
  }
  
  lrtest <- lrtest/iters
  pttest <- pttest/iters
  dptest <- dptest/iters
  
  lrtestind <- lrtestind/iters
  pttestind <- pttestind/iters
  dptestind <- dptestind/iters
  
  return(data.frame(lrtest, pttest, dptest, 
                    lrtestind, pttestind, 
                    dptestind))
}

# function that performs monte-carlo experiment for early hazard difference
EHDExperiment <- function(iters, n=50, censoringrate){
  
  iter <- 0
  lrtest <- 0 # counts number of times log-rank test rejects the null
  pttest <- 0 # counts number of times polya tree prior rejects the null
  dptest <- 0 # counts number of times dirichlet process prior rejects the null
  
  lrtestind <- 0 # counts number of times log-rank test rejects the null
  pttestind <- 0 # counts number of times polya tree prior rejects the null
  dptestind <- 0 # counts number of times dirichlet process prior rejects the null
  
  while(iter < iters){
    print(paste0('Iteration number: ', iter))
    
    group <- c(rep(1, n), rep(2, n)) # group id: treatment group id
    
    survival <- c(rpwexp(n, rate = c(0.75, 3.00, 1.00), intervals = c(0.4, 0.2)), # survival times: treatment group 1
                  rpwexp(n, rate = c(3.00, 0.75, 1.00), intervals = c(0.4, 0.2))) # survival times: treatment group 2
    
    if(censoringrate == 0.00){
      time <- survival
      delta <- rep(1, n)
      data <- data.frame(time, delta, group)
      data <- data[order(data$time),]
    } else{
      if(censoringrate == 0.25){
        b = c(0.621, 0.202)
      } else if(censoringrate == 0.50){
        b = c(1.34, 0.555)
      } else if(censoringrate == 0.75){
        b = c(3.36, 1.921)
      } else{
        print("Invalid Censoring Rate")
        break
      }
      censoring <- c(runif(n, min = 0, max = b[1]), # censoring times: treatment group 1
                     runif(n, min = 0, max = b[2])) # censoring times: treatment group 2
      
      time <- c(NA, n + n) # time: empty vector for event times
      delta <- c(NA, n + n) # delta: survival or censoring time indicator
      
      for(i in 1:(n + n)){
        # time: event times are the first observed event 
        time[i] <- min(survival[i], censoring[i])
        
        # delta: event indicator (1: survival, 0: censoring)
        delta[i] <- ifelse(time[i] == survival[i], 1, 0)
      }
      
      # create the life table for two-sample testing
      data <- data.frame(time, delta, group)
      data <- data[order(data$time),]
    }
    
    # obtain decision from log-rank test and identify if it rejects the null
    lr <- LRTest(data)
    
    # obtain decision from polya-tree test and identify if it rejects the null
    pt <- PTTest(data)
    
    # obtain decision from dirichlet prior test and identify if it rejects the null
    dp <- DPTest(data)
    
    # update counters whenever we get decisive result from dirichlet process prior
    iter <- iter + 1
    
    if(dp == 2){
      lrtestind <- lrtestind + lr
      pttestind <- pttestind + pt
      dptestind <- dptestind + 1
    } else{
      lrtest <- lrtest + lr
      pttest <- pttest + pt
      dptest <- dptest + dp
    }
  }
  
  lrtest <- lrtest/iters
  pttest <- pttest/iters
  dptest <- dptest/iters
  
  lrtestind <- lrtestind/iters
  pttestind <- pttestind/iters
  dptestind <- dptestind/iters
  
  return(data.frame(lrtest, pttest, dptest, 
                    lrtestind, pttestind, 
                    dptestind))
}


# function that performs monte-carlo experiment for late hazard difference
LHDExperiment <- function(iters, n=50, censoringrate){
  
  iter <- 0
  lrtest <- 0 # counts number of times log-rank test rejects the null
  pttest <- 0 # counts number of times polya tree prior rejects the null
  dptest <- 0 # counts number of times dirichlet process prior rejects the null
  
  lrtestind <- 0 # counts number of times log-rank test rejects the null
  pttestind <- 0 # counts number of times polya tree prior rejects the null
  dptestind <- 0 # counts number of times dirichlet process prior rejects the null
  
  
  for(iter in 1:iters){
    group <- c(rep(1, n), rep(2, n)) # group id: treatment group id
    
    survival <- c(rpwexp(n, rate = c(2.00, 0.40), intervals = c(0.50)), # survival times: treatment group 1
                  rpwexp(n, rate = c(2.00, 4.00), intervals = c(0.50))) # survival times: treatment group 2
    
    if(censoringrate == 0.00){
      time <- survival
      delta <- rep(1, n)
      data <- data.frame(time, delta, group)
      data <- data[order(data$time),]
    } else{
      if(censoringrate == 0.25){
        b = c(0.303, 0.303)
      } else if(censoringrate == 0.50){
        b = c(0.92, 0.75)
      } else if(censoringrate == 0.75){
        b = c(3.06, 1.63)
      } else{
        print("Invalid Censoring Rate")
        break
      }
      censoring <- c(runif(n, min = 0, max = b[1]), # censoring times: treatment group 1
                     runif(n, min = 0, max = b[2])) # censoring times: treatment group 2
      
      time <- c(NA, n + n) # time: empty vector for event times
      delta <- c(NA, n + n) # delta: survival or censoring time indicator
      
      for(i in 1:(n + n)){
        # time: event times are the first observed event 
        time[i] <- min(survival[i], censoring[i])
        
        # delta: event indicator (1: survival, 0: censoring)
        delta[i] <- ifelse(time[i] == survival[i], 1, 0)
      }
      
      # create the life table for two-sample testing
      data <- data.frame(time, delta, group)
      data <- data[order(data$time),]
    }
    
    # obtain decision from log-rank test and identify if it rejects the null
    lr <- LRTest(data)
    
    # obtain decision from polya-tree test and identify if it rejects the null
    pt <- PTTest(data)
    
    # obtain decision from dirichlet prior test and identify if it rejects the null
    dp <- DPTest(data)
    
    # update counters whenever we get decisive result from dirichlet process prior
    iter <- iter + 1
    
    if(dp == 2){
      lrtestind <- lrtestind + lr
      pttestind <- pttestind + pt
      dptestind <- dptestind + 1
    } else{
      lrtest <- lrtest + lr
      dptest <- dptest + dp
      pttest <- pttest + pt
    }
  }
  
  lrtest <- lrtest/iters
  pttest <- pttest/iters
  dptest <- dptest/iters
  
  lrtestind <- lrtestind/iters
  pttestind <- pttestind/iters
  dptestind <- dptestind/iters
  

  return(data.frame(lrtest, pttest, dptest, 
                    lrtestind, pttestind, 
                    dptestind))
}



