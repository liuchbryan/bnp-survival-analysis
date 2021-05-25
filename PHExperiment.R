# run code for proportional hazards experiment
require(optparse)

option_list <- list(
  make_option("--rate", default = 0.5,
              help = 'Experimental Censoring Rate')
)

opt <- parse_args(OptionParser(option_list=option_list))

source('BNPExperiments.R')


# run the monte-carlo experiments for each simulated dataset
iters = 100; cr = opt$rate
results <- PHExperiment(iters = iters, n = 50, censoringrate = cr)

cr <- as.integer(cr*100)

filename <- paste0('PHExperiment', cr, '.rds')
saveRDS(results,file=filename)