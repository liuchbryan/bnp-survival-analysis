require(optparse)

option_list <- list(
  make_option("--rate", default = 0.5,
              help = 'Experimental Censoring Rate')
)

opt <- parse_args(OptionParser(option_list=option_list))

source('BNPExperiments.R')


# run the monte-carlo experiments for each simulated dataset
iters = 100; cr = opt$rate
results <- SHExperiment(iters = iters, n = 25, censoringrate = cr)
# PHExperiment(iters = iters, n = 25, censoringrate = cr)
# EHDExperiment(iters = iters, n = 25, censoringrate = cr)
# LHDExperiment(iters = iters, n = 25, censoringrate = cr)

cr <- gsub('\\.', '_', cr)

filename <- paste0('Experiment_SH', cr, '.rds')
saveRDS(results,file=filename)