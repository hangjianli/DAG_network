# load packages----------------------------------------------------------------
rm(list = ls())
setwd("~/Documents/research/dag_network/")
source("R/newalgo/loadpackages.R")
source("R/newalgo/helper_funcs.R")
source("R/newalgo/gen_params_funcs.R")
source("R/newalgo/functions.R")
source("R/newalgo/flipflop.R")
options(scipen = 7)

# generate specs ----------------------------------------------------------
args <- args_for_parameter()  
dir.create(path = paste0('output/', args$setting))
saveRDS(args, file = paste0('output/', args$setting, "/args.rds"))
estimands <- generate_parameters(
  args = args, 
  # bname = 'hepar2',
  # btype = 'discrete',
  # theta_name = 'facebook',
  sbm_nblocks = 5,
  cluster_sizes = c(20, 15, 5, 5, 10, 10, 20, 20, 15, 25, 15, 5, 15, 20),
  sbm_probs = c(0.01,0.5),
  seed = 3)

image(estimands$theta, axes=T)
image(as(estimands$b, class(estimands$theta)))
estimands$b[estimands$b!=0] %>% as.numeric() %>% hist(breaks=20)
estimands$theta[estimands$theta!=0] %>% as.numeric() %>% hist(breaks=20)
# main_lam_6$thetahat[main_lam_6$thetahat!=0] %>% as.numeric %>% hist(breaks=20)
saveRDS(estimands, file = paste0('output/', args$setting, "/estimands.rds"))
estimands$s0
# X_ <- sim_X(
#   vers = 1,
#   n = args$n,
#   p = estimands$realp,
#   omg.sq = estimands$omg.sq,
#   sig = estimands$sig,
#   b = estimands$b
# )
# run simulation ----------------------------------------------------------
# ordered 
setwd("~/Documents/research/dag_network")
simID <- args$setting
sim_newalgo_ordered(
  args,
  estimands, 
  start_sim=1,
  end_sim=1, 
  lamLen=10,
  lambda2 = 1,
  lambda1_max_div = 10,
  lambda1_max_div2 = 250)

# as(main_lam_6$thetahat[1:10,1:10], 'dgCMatrix')
# image(estimands$theta, axes=T)
eigen(estimands$theta %>% round(digits = 2))$value

# simID = '121006'
# args <- readRDS(paste0("~/Documents/research/dag_network/output/", simID,"/args.rds"))
# estimands <- readRDS(paste0("~/Documents/research/dag_network/output/", simID, "/estimands.rds"))
# process_output_ordered_tmp(
#   simID=args$setting,
#   estimands,
#   args,
#   start=5,
#   num_sim=10,
#   thr=0.1
# )
# get_average_shd_ordered(simID = args$setting, nsim = 10)
# 

simID = args$setting
process_output_ordered(simID = args$setting, estimands = estimands, 
                       args = args, start = 1, num_sim = 1,  thr = 0.1, ff_flag = T)
get_all_shd_ordered(
  simID = args$setting, 
  estimands,
  start = 1, 
  num_sim = 1,
  ff_flag = T
)
get_average_shd_ordered(simID = args$setting, nsim = 1, ff_flag = T)

setwd("~/Documents/research/dag_network")
get_roc_plots_noff(simID = args$setting, 
                   ymin = 20, 
                   xmax = 150,
                   thrLen = 50,
                   start = 1, 
                   nsim = 1,
                   ff_flag = T)

# shd_average <- readRDS("~/Documents/research/dag_network/output/101/shd_average.rds")
# shd_average$shdXbaseline
# shd_average$shdXff
args <- readRDS(paste0("~/Documents/research/dag_network/output/", simID,"/args.rds"))
estimands <- readRDS(paste0("~/Documents/research/dag_network/output/", simID, "/estimands.rds"))
# get_roc_data_ordered(simID = simID)
get_all_shd_ordered(simID = simID, estimands, args$num_sim)
get_average_shd_ordered(simID = simID, nsim = as.numeric(args$num_sim))

simID = '8003'
shd_average <- readRDS(paste0("~/Documents/research/dag_network/output/", simID,"/shd_average.rds"))
xtable(t(shd_average), digits = c(0,1,1,1,3,3,1,3,3,1))
shd_average

