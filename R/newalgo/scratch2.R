setwd("~/Documents/research/dag_network")
# args <- readRDS("~/Documents/research/dag_network/output/121006/args.rds")
# estimands <- readRDS("~/Documents/research/dag_network/output/121006/estimands.rds")
nsim = args$num_sim
# nsim = 1
for(sim in 2:2){
  # dir.create(path = paste0("output/",args$setting, "/", args$setting, "--", sim))
  setwd(paste0("output/",args$setting, "/", args$setting, "--", sim))
  con <- file("test2.log")
  sink(con, append=TRUE)
  sink(con, append=TRUE, type="message")
  cat("===============\n")
  cat(paste0('starting simulation ', sim, '\n'))
  X_ <- readRDS(paste0("~/Documents/research/dag_network/output/",args$setting, "/", 
                       args$setting, "--", sim, '/X.rds'))
  X <- X_$X
  networkDAG_sol_path2(
    X = X,
    block_size = args$block_size,
    zeropos_list = estimands$zeropos_list,
    lambda_len = 1,
    lambdaff_max = 50,
    lambdaff_min = 10,
    maxIter = 4,
    lambda2 = 100
  )
  sink() 
  sink(type="message")
  setwd("~/Documents/research/dag_network")
}
simID = args$setting
process_output_ordered(simID = simID, estimands = estimands, thr = 0.1)
get_all_shd_ordered(simID = simID, estimands, args$num_sim)
get_average_shd_ordered(simID = simID, nsim = as.numeric(args$num_sim))

