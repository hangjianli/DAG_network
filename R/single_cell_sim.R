sig_genes_log_val <- readRDS("data/sig_genes_log_val_full.rds")
goodgenes_names <- unlist(read_excel("data/single_cell_data/single_cell.xlsx"))
goodgene <- sig_genes_log_val[rownames(sig_genes_log_val) %in% goodgenes_names, ]
