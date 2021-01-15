setwd("~/Documents/research/dag_network/")
library(dplyr)
targetgene <- readRDS("data/single_cell_data/sig_genes_log_val.rds")
sig_genes_log_val <- readRDS("data/sig_genes_log_val_full.rds")
idx_goodgene <- apply(sig_genes_log_val, 1, function(x) sd(x)/abs(mean(x)) > 0.25) %>% which()
goodgene <- sig_genes_log_val[idx_goodgene,]
goodgene %>% dim()
set.seed(12)
othergenes <- sample_n(as.data.frame(goodgene[1:7000,]), 2000) %>% as.matrix() %>% t()



# default clustering ------------------------------------------------------
sc_block_idx_full <- readRDS("data/single_cell_data/single_cell_block_idx_full.rds")
for(i in 1:length(sc_block_idx_full)){
  names(sc_block_idx_full[[i]]) <- colnames(targetgene)[sc_block_idx_full[[i]]]
}
# Xp <- readRDS(file = "data/single_cell_data/sig_genes_log_val.rds")
Xp <- t(targetgene)
Xp %>% dim()
# randomly sample 20 cells from each cell type and merge the indices
res <- sample_sc_data(
  full_log_vals = Xp, 
  full_idx = sc_block_idx_full,
  size = c(10,10,20,10,15,15,15),
  seed = 4
)
# res <- Xp
res$subsetXp %>% dim()
saveRDS(res, "sc_subsampled.rds")


# clustering --------------------------------------------------------------
d <- dist(scale(othergenes), method = "euclidean")
hc1 <- hclust(d, method = "complete")
# plot(hc1, cex = 0.6, hang = -1)
sub_grp <- cutree(hc1, h=75)
sub_grp %>% table()
sub_grp_subset <- sub_grp[!(sub_grp %in%  which(table(sub_grp) ==1))]
sub_grp_subset %>% table()
sub_grp_subset %>% table() %>% max()
sub_grp_remove <- sub_grp[(sub_grp %in%  which(table(sub_grp) == 1))]
sub_grp_subset %>% length()
# 
# test <- substr(names(sub_grp_subset), start = 0, stop = 3)
# cell_names <- sub(pattern = "_$", "", x = test)[-1] %>% unique()
# for(n in cell_names){
#   i_H1 <- grepl(paste0("^", n), names(sub_grp))
#   H1 <- sub_grp[i_H1]
#   cat("[INFO]: ", n)
#   H1 %>% table() %>% print()
# }


# targetgene <- t(targetgene)
targetgene %>% dim()
res <- reorder_data(sub_grp_subset = sub_grp_subset, targetgene = targetgene)
res$df %>% dim()
res$block_idx[[1]]
res$block_idx %>% length()



# run sim -----------------------------------------------------------------

dir.create(path = "output/single_cell13")
setwd(dir = "output/single_cell13")

networkDAG_sol_path(
  # X = res$subsetXp,
  X = res$df,
  block_size=20, 
  zeropos_list = NULL,
  block_idx = res$block_idx,
  lambda_len = 10,
  maxIter = 100
)
df <- res$df

Xdecor_res <- get_Xdecor(df)

GES_sol(df, decor = F)
GES_sol(Xdecor_res$X_decor, decor = T)
# GES_sol(Xdecor_res$X_decor_1iter, decor = T)
pc_sol(df, decor = F)
pc_sol(Xdecor_res$X_decor, decor = T)
sparsebn_sol(df, decor = F)
sparsebn_sol(Xdecor_res$X_decor, decor = T)

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