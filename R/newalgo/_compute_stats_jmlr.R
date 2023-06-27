rm(list = ls())
setwd("~/Documents/research/dag_network/")
source("R/newalgo/loadpackages.R")
source("R/newalgo/helper_funcs.R")
source("R/newalgo/gen_params_funcs.R")
source("R/newalgo/functions.R")
source("R/newalgo/flipflop.R")
options(scipen = 7)

start_sim=1
end_sim=10

setwd("~/Documents/research/dag_network")
simID <- '003'
args <- readRDS(paste0("~/Documents/research/dag_network/output/", simID,"/args.rds"))
estimands <- readRDS(paste0("~/Documents/research/dag_network/output/", simID, "/estimands.rds"))
start_time = Sys.time()
# 
# for(sim in start_sim:end_sim){
#   setwd(paste0("output/",args$setting, "/", args$setting, "--", sim))
#   cat("------------------\n")
#   cat(paste0("Start sim ", sim, " at ", round(difftime(Sys.time(), start_time, units = "mins")[[1]], 2), " minutes \n"))
#   # con <- file("test_jmlr.log")
#   # sink(con, append=TRUE)
#   # sink(con, append=TRUE, type="message")
#   X_ = readRDS("X.rds")
#   X <- X_$X
#   res = bnlearn::pc.stable(
#     x = as.data.frame(X),
#     blacklist=ordering2blacklist(dimnames(X)[[2]]),
#     alpha = 0.1,
#   )
#   cat(paste0("Encode structural constraints", " at ", round(difftime(Sys.time(), start_time, units = "mins")[[1]], 2), " minutes. \n"))
#   fitted = bn.fit(x = res, data = as.data.frame(X), method = 'mle-g')
#   saveRDS(fitted, file = paste0("pc_ordered.rds"))
#   # estimate MLE for PC -----------------------------------------------------
#   n = dim(X)[1]
#   p = dim(X)[2]
#   pc_res <- readRDS(file = 'pc_ordered.rds')
#   pc_adj_b = matrix(0,p,p)
#   for(i in 2:p){
#     if(length(pc_res[[i]]$coefficients) > 1){
#       pc_adj_b[as.integer(names(pc_res[[i]]$coefficients)[-1]), 
#                i] = pc_res[[i]]$coefficients[-1]
#     }
#   }
#   mle_result <- dag_mle_estimation(
#     X = X,
#     Bhat = pc_adj_b,
#     Lhat = diag(n)
#   )
#   saveRDS(mle_result, file = paste0("pc_MLE_result.rds"))
#   # BIC_result <- BIC_dag(
#   #   X = X,
#   #   block_idx = block_idx,
#   #   bmle = mle_result$Bmle,
#   #   omgmle = mle_result$omgmlesq,
#   #   theta = diag(n)
#   # )
#   # saveRDS(BIC_1iter_result,  paste0('BIC_1iter_result_', k, '.rds'))
#   # sink()
#   # sink(type="message")
#   cat(paste0("Experiment ", sim, " completed", " at ", round(difftime(Sys.time(), start_time, units = "mins")[[1]], 2), " minutes. \n"))
#   setwd("~/Documents/research/dag_network")
# }
# 
# 
setwd("~/Documents/research/dag_network")



simID = args$setting
process_output_ordered(
  simID = args$setting,
  estimands = estimands, 
  args = args, 
  start = 1, 
  num_sim = 10,  
  threshholds = 0.1, 
  ff_flag = T,
  pc_flag = T,
  ges_flag = T
  # threshholds=seq(0, 0.5,length.out = 50)
)

# 

get_all_shd_ordered(
  stats_file_name = 'SHDstats_JMLR.rds',
  simID = args$setting, 
  estimands,
  start = 1, 
  num_sim = 10,
  ff_flag = T, 
  pc_flag = T, 
  ges_flag = T
)


get_average_shd_ordered(simID = args$setting, nsim = 10, ff_flag = T, pc_flag = T, ges_flag = T)


