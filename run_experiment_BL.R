require(optparse)

option_list <- list(
  make_option("--rate", default = 0.5,
              help = 'Experimental Censoring Rate')
)

opt <- parse_args(OptionParser(option_list=option_list))

source('BNPExperiments.R')


# run the monte-carlo experiments for each simulated dataset
iters = 100; cr = opt$rate
#Andrea
#results <- SHExperiment(iters = iters, n = 25, censoringrate = cr)
# Ben
# results <- PHExperiment(iters = iters, n = 25, censoringrate = cr)
# Bryan 
results <-  EHDExperiment(iters = iters, n = 25, censoringrate = cr)
# Jose
# results <- LHDExperiment(iters = iters, n = 25, censoringrate = cr)

cr <- gsub('\\.', '-', cr)

filename <- paste0('Experiment_EHD_', n, '_', cr, '.rds')
saveRDS(results,file=filename)