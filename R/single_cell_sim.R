library(dplyr)
targetgene <- readRDS("data/single_cell_data/sig_genes_log_val.rds")
sig_genes_log_val <- readRDS("data/sig_genes_log_val_full.rds")
idx_goodgene <- apply(sig_genes_log_val, 1, function(x) sd(x)/abs(mean(x)) > 0.25) %>% which()
goodgene <- sig_genes_log_val[idx_goodgene,]
goodgene %>% dim()
set.seed(10)
othergenes <- sample_n(as.data.frame(goodgene[1:7000,]), 3000) %>% as.matrix() %>% t()





# clustering --------------------------------------------------------------
d <- dist(scale(othergenes), method = "euclidean")
hc1 <- hclust(d, method = "complete")
# plot(hc1, cex = 0.6, hang = -1)
sub_grp <- cutree(hc1, h=65)
sub_grp %>% table()
sub_grp_subset <- sub_grp[!(sub_grp %in%  which(table(sub_grp) == 1))]
sub_grp_subset %>% table()
sub_grp_remove <- sub_grp[(sub_grp %in%  which(table(sub_grp) == 1))]
sub_grp_remove

# 
# test <- substr(names(sub_grp_subset), start = 0, stop = 3)
# cell_names <- sub(pattern = "_$", "", x = test)[-1] %>% unique()
# for(n in cell_names){
#   i_H1 <- grepl(paste0("^", n), names(sub_grp))
#   H1 <- sub_grp[i_H1]
#   cat("[INFO]: ", n)
#   H1 %>% table() %>% print()
# }

if(any(table(sub_grp_subset) == 1)){
  warnings('Too many blocks!')
}
group_idx = unique(sub_grp_subset)
block_idx = vector(mode = 'list', length = length(group_idx))
for(i in 1:length(block_idx)){
  block_idx[[i]] = which(sub_grp_subset == group_idx[i])
}



dir.create(path = "output/single_cell11")
setwd(dir = "output/single_cell11")

targetgene <- readRDS("~/Documents/research/dag_network/data/single_cell_data/sig_genes_log_val.rds")
targetgene <- targetgene[, !(colnames(targetgene) %in% names(sub_grp_remove))]
Xp <- t(targetgene)
Xp %>% dim()
networkDAG_sol_path(
  # X = res$subsetXp, 
  X = Xp,
  block_size=20, 
  zeropos_list = NULL,
  block_idx = block_idx,
  lambda_len = 10,
  maxIter = 100
)
