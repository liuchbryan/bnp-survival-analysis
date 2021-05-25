# run code for late difference experiment
require(optparse)

option_list <- list(
  make_option("--rate", default = 0.5,
              help = 'Experimental Censoring Rate')
)

opt <- parse_args(OptionParser(option_list=option_list))

source('BNPExperiments.R')

# run the monte-carlo experiments for each simulated dataset
iters = 100; cr = opt$rate
results <- LHDExperiment(iters = iters, n = 50, censoringrate = cr)

cr <- as.integer(cr*100)

filename <- paste0('LDExperiment', cr, '.rds')
saveRDS(results,file=filename)