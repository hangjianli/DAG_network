# real data ---------------------------------------------------------------
library(readr)
library(dplyr)
library(readxl)
# single_cell <- read_csv("data/single_cell_data/GSE75748_sc_cell_imputationscimpute_count.csv")
# single_cell %>% dim()
# names(single_cell)[1] <- 'geneID'
# saveRDS(single_cell, "data/single_cell.rds")

# single_cell <- readRDS("data/single_cell_data/single_cell.rds")
# sig_genes <- single_cell
# %>% filter(geneID %in% gene_names$GeneID)
# sig_genes %>% dim()
# count freq of each cell type
# test <- substr(names(sig_genes), start = 0, stop = 3)
# sub(pattern = "_$", "", x = test)[-1] %>% table()
# log transform of gene count data
# sig_genes[, -1] <- log10(sig_genes[, -1] + 1.01)
# sig_genes_log_val <- sig_genes[, -1] %>% as.matrix()
# dimnames(sig_genes_log_val)[[1]] <- sig_genes$geneID
# saveRDS(sig_genes_log_val, "data/sig_genes_log_val_full.rds")


# -------------------------------------------------------------------------
# goodgene consists of 51 target genes and 1018 cells
# othergenes consists of 3000 genes and 1018 cells used to estimate blocks
# -------------------------------------------------------------------------

sig_genes_log_val <- readRDS("data/sig_genes_log_val_full.rds")
goodgenes_names <- unlist(read_excel("data/single_cell_data/single_cell.xlsx"))
# find set of valid genes by coefficient of variation ( sd / |mean| > 0.25)
idx_goodgene <- apply(sig_genes_log_val, 1, function(x) sd(x)/abs(mean(x)) > 0.25) %>% which()
idx_goodgene %>% length()
goodgene <- sig_genes_log_val[idx_goodgene,]
goodgene[1:5,1:5]
set.seed(1)
othergenes <- sample_n(as.data.frame(goodgene[1:7000,]), 3000) %>% as.matrix() %>% t()
othergenes %>% dim()
# select 51 good genes ----------------------------------------------------
goodgene <- sig_genes_log_val[rownames(sig_genes_log_val) %in% goodgenes_names, ]
goodgene %>% str()

# hierarchical clustering here---------------------------------------------
res$subsetXp %>% dim()
d <- dist(scale(res$subsetXp), method = "euclidean")
hc1 <- hclust(d, method = "complete" )
plot(hc1, cex = 0.6, hang = -1)
sub_grp <- cutree(hc1, k = 10)
table(sub_grp)




res$cellnames







# reference ---------------------------------------------------------------

# compute sample covariance for each cell type ----------------------------
i_H1 <- grepl("^H1", row.names(othergenes))
H1 <- othergenes[i_H1,]
test <- cov(t(H1))
test %>% dim()

test %>% range()

solve(test)[1:10,1:10]

i_DEC <- grepl("^DEC", row.names(othergenes))
i_EC <- grepl("^EC", row.names(othergenes))
i_H9 <- grepl("^H9", row.names(othergenes))
i_HFF <- grepl("^HFF", row.names(othergenes))
i_NPC <- grepl("^NPC", row.names(othergenes))
i_TB <- grepl("^TB", row.names(othergenes))


