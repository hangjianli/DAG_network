# load packages----------------------------------------------------------------
rm(list = ls())
setwd("~/Documents/research/dag_network/")
source("R/loadpackages.R")
source("R/helper_funcs.R")
source("R/gen_params_funcs.R")
source("R/newalgo/functions.R")
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
setwd("~/Documents/research/dag_network")
sim_newalgo_ordered(args, estimands, start_sim=1, end_sim=args$num_sim, lamLen=15)
process_output_ordered(simID='003', thr = 0.1)
get_all_shd_ordered(simID = '003', estimands = estimands, num_sim = args$num_sim)
get_average_shd(simID = '003', nsim = as.numeric(args$num_sim))





sim_newalgo_unordered <- function(
  args, 
  estimands, 
  start_sim=1, 
  end_sim=args$num_sim, 
  lamLen=15){
  
  
}