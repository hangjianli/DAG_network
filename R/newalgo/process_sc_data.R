# real data ---------------------------------------------------------------
library(readr)
# single_cell <- read_csv("data/single_cell_data/GSE75748_sc_cell_imputationscimpute_count.csv")
# single_cell %>% dim()
# names(single_cell)[1] <- 'geneID'
# saveRDS(single_cell, "data/single_cell.rds")

# 
single_cell <- readRDS("data/single_cell_data/single_cell.rds")
sig_genes <- single_cell 
# %>% filter(geneID %in% gene_names$GeneID)
sig_genes %>% head()
# test <- substr(names(sig_genes), start = 0, stop = 3)
# sub(pattern = "_$", "", x = test)[-1] %>% table() 
sig_genes[, -1] <- log10(sig_genes[, -1] + 1.01)
sig_genes_log_val <- sig_genes[, -1] %>% as.matrix()
dimnames(sig_genes_log_val)[[1]] <- sig_genes$geneID
saveRDS(sig_genes_log_val, "data/sig_genes_log_val_full.rds")

# sig_genes_log_val %>% dim()
library(dplyr)

# find set of valid genes by coefficient of variation

idx_goodgene <- apply(sig_genes_log_val, 1, function(x) sd(x)/abs(mean(x)) > 0.2) %>% which()

goodgene <- sig_genes_log_val[idx_goodgene,]
goodgene_3k <- sample_n(as.data.frame(goodgene), 3000) %>% as.matrix()
goodgene_3k %>% str()
# hierarchical clustering here
library(hclust)
