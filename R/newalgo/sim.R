# load packages----------------------------------------------------------------
rm(list = ls())
setwd("~/Documents/research/dag_network/")
source("R/loadpackages.R")
source("R/helper_funcs.R")
source("R/gen_params_funcs.R")
source("R/newalgo/functions.R")
source("R/newalgo/flipflop.R")


# generate specs ----------------------------------------------------------
args <- args_for_parameter()  
dir.create(path = paste0('output/', args$setting))
saveRDS(args, file = paste0('output/', args$setting, "/args.rds"))
estimands <- generate_parameters(args = args, seed = 1, )
image(estimands$theta)
image(as(estimands$b, class(estimands$theta)))
estimands$b[estimands$b!=0] %>% as.numeric() %>% hist(breaks=20)
saveRDS(estimands, file = paste0('output/', args$setting, "/estimands.rds"))
# run simulation ----------------------------------------------------------
# ordered 
setwd("~/Documents/research/dag_network")
simID <- args$setting
sim_newalgo_ordered(
  args,
  estimands, 
  start_sim=1,
  end_sim=args$num_sim, 
  lamLen=10,
  lambda1_max_div = 100)
process_output_ordered(simID = simID, estimands = estimands, thr = 0.1)
get_all_shd_ordered(simID = simID, estimands, args$num_sim)
get_average_shd_ordered(simID = simID, nsim = as.numeric(args$num_sim))

simID = '121005'


shd_average <- readRDS("~/Documents/research/dag_network/output/121006/shd_average.rds")
shd_average$shdXbaseline
shd_average$shdXff


args <- readRDS(paste0("~/Documents/research/dag_network/output/", simID,"/args.rds"))
estimands <- readRDS(paste0("~/Documents/research/dag_network/output/", simID, "/estimands.rds"))
# get_roc_data_ordered(simID = simID)
get_all_shd_ordered(simID = simID, estimands, args$num_sim)
get_average_shd_ordered(simID = simID, nsim = as.numeric(args$num_sim))
get_roc_plots(simID = simID, ymin = 20, xmax = 150)


# unordered 
setwd("~/Documents/research/dag_network")
simID <- args$setting
sim_newalgo_unordered(args, estimands, start_sim=4, end_sim=args$num_sim, lamLen=10, lambda1_max_div = 100)
process_output_unordered(simID = simID, thr = 0.1, nsim = args$num_sim)
get_average_shd_unordered(simID = simID, nsim = as.numeric(args$num_sim))


simIDarr = c(
  '00011','00022', '00033', '00055', '00044',
  '321001', '321002', '321003', '321004', '321005'
)

shd_data <- get_shd_ji_diff(simIDarr)
plot_boxplot_shd_ja(shd_data = shd_data)


# library(xtable)
# xtable(shd_average, digits = c(0,0,0))


