rm(list = ls())
setwd("~/Documents/research/dag_network/")
source("R/newalgo/loadpackages.R")
source("R/newalgo/helper_funcs.R")
source("R/newalgo/gen_params_funcs.R")
source("R/newalgo/functions.R")
source("R/newalgo/flipflop.R")
options(scipen = 7)


# n < p -------------------------------------------------------------------

# Rerun PC and GES with given order, assuming observations are independent.
# 101, 002, 003, 007

start_sim=1
end_sim=10

setwd("~/Documents/research/dag_network")
simID <- '101'
args <- readRDS(paste0("~/Documents/research/dag_network/output/", simID,"/args.rds"))
estimands <- readRDS(paste0("~/Documents/research/dag_network/output/", simID, "/estimands.rds"))
start_time = Sys.time()

get_roc_data_ordered(simID = simID, thrub = 1.2)

get_roc_plots_noff(
  simID = args$setting, 
  ymin = 20, 
  xmax = 150,
  thrLen = 20,
  start = 1, 
  nsim = 10,
  ff_flag = T
)

