library(dplyr)

setwd("~/Documents/research/dag_network")
simID <- '121002'
args <- readRDS(paste0("~/Documents/research/dag_network/output/", simID,"/args.rds"))
estimands <- readRDS(paste0("~/Documents/research/dag_network/output/", simID, "/estimands.rds"))
# get_roc_data_ordered(simID = simID)


# get_all_shd_ordered(simID = simID, estimands, num_sim=args$num_sim, ff_flag = T)
get_average_shd_ordered(simID = simID, nsim = as.numeric(args$num_sim), ff_flag = T)
shd_average_jmlr <- readRDS(paste0("~/Documents/research/dag_network/output/", simID,"/shd_average_jmlr.rds"))
shd_se <- readRDS(paste0("~/Documents/research/dag_network/output/", simID,"/shd_se_jmlr.rds"))
xtable(t(shd_average_jmlr), digits = c(0,1,1,1,3,3,1,3,3,1))
# shd_average
# shd_se

xtable(t(shd_se), digits = c(0,1,1,1,3,3,1,3,3,1))

# aggregate(. ~ ID, DF, function(x) c(mean = mean(x), sd = sd(x)))
# shd_average_old_results <- readRDS(paste0("~/Documents/research/dag_network/output/", simID,"/shd_average.rds"))
# xtable(t(shd_average_old_results), digits = c(0,1,1,1,3,3,1,3,3,1))

# 
# 
# res = c()
# for (sim in 1:10){
#   allshd <- as.data.frame(
#     readRDS(paste0("output/", simID,  '/', simID, "--", sim, "/SHDclose_jmlr.rds"))
#   )
#   res = c(res, allshd['pnum', 'shdXmain'])
# }


# unordered ---------------------------------------------------------------


setwd("~/Documents/research/dag_network")
simID <- '603'
args <- readRDS(paste0("~/Documents/research/dag_network/output/", simID,"/args.rds"))
estimands <- readRDS(paste0("~/Documents/research/dag_network/output/", simID, "/estimands.rds"))
# get_roc_data_ordered(simID = simID)


# get_all_shd_ordered(simID = simID, estimands, num_sim=args$num_sim, ff_flag = T)
get_average_shd_unordered(start = 1, simID = simID, nsim = as.numeric(args$num_sim))


# manual check ------------------------------------------------------------
test = numeric(10)
for (sim in 1:10){
  allshd <- as.data.frame(
    readRDS(paste0("output/", simID,  '/', simID, "--", sim, "/SHDclose.rds"))
  )
  print(allshd['pnum', 'shdXmain1iter'])
  test[sim] = allshd['pnum', 'shdXmain1iter']
  # allshd <- allshd[, names(total)]
  # allshd_list[[sim]] <- allshd
  # total <- total + allshd
}
allshd_list = list()
for (sim in 1:nsim){
  allshd <- as.data.frame(
    readRDS(paste0("output/", simID,  '/', simID, "--", sim, "/SHDclose.rds"))
  )
  allshd_list[[sim]] <- allshd
}
