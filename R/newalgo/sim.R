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
sim_newalgo_ordered(args, estimands, start_sim=1, end_sim=args$num_sim, lamLen=10)
process_output_ordered(simID = simID, estimands = estimands, thr = 0.1)
get_all_shd_ordered(simID = simID, estimands, args$num_sim)
get_average_shd_ordered(simID = simID, nsim = as.numeric(args$num_sim))

# unordered 
setwd("~/Documents/research/dag_network")
simID <- args$setting
sim_newalgo_unordered(args, estimands, start_sim=1, end_sim=10, lamLen=10)
process_output_unordered(simID = simID, thr = 0.1, nsim = args$num_sim)
get_average_shd_unordered(simID = simID, nsim = as.numeric(args$num_sim))

# single cell data --------------------------------------------------------




# hierarchical clustering here---------------------------------------------
num_blocks = 5
res$subsetXp %>% dim()

# outsidedata <- othergenes[rownames(othergenes) %in% rownames(res$subsetXp), ]
outsidedata <- othergenes[rownames(othergenes) %in% rownames(res), ]

d <- dist(scale(outsidedata), method = "euclidean")
hc1 <- hclust(d, method = "complete")
# plot(hc1, cex = 0.6, hang = -1)
sub_grp <- cutree(hc1, k = num_blocks)
table(sub_grp)
if(any(table(sub_grp) == 1)){
  warnings('Too many blocks!')
}
block_idx = vector(mode = 'list', length = num_blocks)
for(i in 1:num_blocks){
  block_idx[[i]] = which(sub_grp == i)
}


# 
# sum(sapply(block_idx, length) <= 2)
# 
# block_idx2 = vector(mode = 'list')
# 
# i = 1
# zeropos_list <- matrix(0, dimnames = list(NULL, c('row', 'col')), ncol = 2)
# newblock = c()
# for(j in 1:k){
#   if(length(block_idx[[j]]) > 2){
#     block_idx2[[i]] = block_idx[[j]]
#     i <- i + 1
#   }else{
#     newblock <- c(newblock, block_idx[[j]])
#     zeropos_list <- rbind(zeropos_list, )
#   }
# }


# test_mat <- matrix(1:16, 4, 4)
# dimnames(test_mat) <- list(1:4, 1:4)
networkDAG_sol_path(
  X = res$subsetXp,
  # X = Xp,
  block_size=20, 
  zeropos_list = NULL,
  block_idx = res$block_idx,
  lambda_len = 10,
  maxIter = 100
)


# 
# X: num cells x num genes
# X = [X1, X2]
# 
# X2 -> theta_hat
# X1 -> cross validation 
# 
# likelihood1: CV(X1) without decor using PC, GES, ... theta_hat
# likelihood2: CV(X1) with decor using PC, GES, ... theta_hat
# 
# in first case:
#   X1 = X11 -> B_hat1,  X12 + B_hat1 + theta_hat -> likelihood1
# 
# in second case:
#   X1 = decor(X11) -> B_hat2,  X12 + B_hat2 + theta_hat -> likelihood2
# 
# 
