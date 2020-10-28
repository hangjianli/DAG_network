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
# ordered 
setwd("~/Documents/research/dag_network")
simID <- args$setting
sim_newalgo_ordered(args, estimands, start_sim=1, end_sim=args$num_sim, lamLen=15)
process_output_ordered(simID = simID, thr = 0.1)
get_all_shd_ordered(simID = simID, estimands, args$num_sim)
get_average_shd_ordered(simID = simID, nsim = as.numeric(args$num_sim))

# unordered 
setwd("~/Documents/research/dag_network")
simID <- args$setting
sim_newalgo_unordered(args, estimands, start_sim=1, end_sim=10, lamLen=15)
process_output_unordered(simID = simID, thr = 0.1)
get_average_shd_unordered(simID = simID, nsim = as.numeric(args$num_sim))
