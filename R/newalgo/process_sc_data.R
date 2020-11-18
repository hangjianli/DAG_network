


# real data ---------------------------------------------------------------

single_cell <- read_csv("data/GSE75748_sc_cell_imputationscimpute_count.csv")
single_cell %>% dim()
names(single_cell)[1] <- 'geneID'
# saveRDS(single_cell, "data/single_cell.rds")
single_cell[1:5,1:5]

library(readxl)
gene_names <- read_excel("data/single_cell.xlsx")
gene_names

"EOMES" %in% single_cell$geneID

# 
sig_genes <- single_cell %>% filter(geneID %in% gene_names$GeneID)

sig_genes %>% head()

test <- substr(names(sig_genes), start = 0, stop = 3)
sub(pattern = "_$", "", x = test)[-1] %>% table() 
sig_genes[, -1] <- log10(sig_genes[, -1] + 1.01)


sig_genes_log_val <- sig_genes_log[, -1] %>% as.matrix()

dimnames(sig_genes_log_val)[[1]] <- sig_genes$geneID

saveRDS(sig_genes_log_val, "data/sig_genes_log_val.rds")
