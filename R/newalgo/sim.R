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

setwd("~/Documents/research/dag_network/")
sc_block_idx_full <- readRDS("data/single_cell_data/single_cell_block_idx_full.rds")
# Xp <- readRDS(file = "data/single_cell_data/sig_genes_log_val.rds")
dir.create(path = "output/single_cell03")
setwd(dir = "output/single_cell01")
Xp <- t(goodgene)
Xp %>% dim()
set.seed(2)
# randomly sample 20 cells from each cell type and merge the indices
res <- sample_sc_data(
  full_log_vals = Xp, 
  full_idx = sc_block_idx_full,
  size = c(20,20,20,10,15,15,15),
  seed = 4
)
res <- Xp
res$subsetXp %>% dim()
saveRDS(res, "sc_subsampled.rds")



# hierarchical clustering here---------------------------------------------
num_blocks = 5
# res$subsetXp %>% dim()

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

# get the block index -----------------------------------------------------
# block_idx <- get_cell_block_idx(cellnames = res$cellnames)


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





# combine small blocks -----------------------------------------------------



# correct_order <- match(rownames(res$subsetXp), names(unlist(block_idx)))


# test_mat <- matrix(1:16, 4, 4)
# dimnames(test_mat) <- list(1:4, 1:4)
networkDAG_sol_path(
  # X = res$subsetXp, 
  X = Xp,
  block_size=20, 
  zeropos_list = NULL,
  block_idx = block_idx,
  lambda_len = 10,
  maxIter = 100
)
Xdecor_res <- get_Xdecor(res$subsetXp)

GES_sol(res$subsetXp, decor = F)
GES_sol(Xdecor_res$X_decor, decor = T)
# GES_sol(Xdecor_res$X_decor_1iter, decor = T)
pc_sol(Xdecor_res$X_decor, decor = T)
pc_sol(res$subsetXp, decor = F)
sparsebn_sol(Xdecor_res$X_decor, decor = T)
sparsebn_sol(res$subsetXp, decor = F)
fgesdag <- readRDS("adjmat_fges_CPDAG_decor.rds")
fgesdag_original <- readRDS("adjmat_fges_CPDAG.rds")
pcdag <- readRDS("adjmat_pc_CPDAG_decor.rds")
pcdag_original <- readRDS("adjmat_pc_CPDAG.rds")
sbndag <- readRDS("adjmat_sparsebn_CPDAG_decor.rds")
sbndag_original <- readRDS("adjmat_sparsebn_CPDAG.rds")

plot_cpdag(fgesdag_original)
plot_cpdag(fgesdag)
plot_cpdag(pcdag_original)
plot_cpdag(pcdag)
plot_cpdag(sbndag_original)
plot_cpdag(sbndag,rescale = T)



bic_score <- readRDS('BICscores_main.rds')
best_res <- readRDS(paste0('main_lam_', best_bic, '.rds'))
bstar_adj_cpdag <- bnstruct::dag.to.cpdag(1*(best_res$bhat != 0)) 
which(best_res$bhat != 0, arr.ind = T)
plot_cpdag(best_res$bhat)
plot_cpdag(best_res$bhat_1iter)


which(best_res$bhat[, 'POU5F1'] != 0)

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
