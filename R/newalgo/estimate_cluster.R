setwd("~/Documents/research/dag_network/")

targetgene <- readRDS("data/single_cell_data/sig_genes_log_val.rds")
targetgene <- t(targetgene)
targetgene %>% dim() # 51 target genes
goodgene <- readRDS('data/single_cell_data/allgoodgenes.rds')
goodgene %>% dim() # all the genes that passed VOC test
goodgene_exclude_target <- goodgene[!(rownames(goodgene) %in% colnames(targetgene)), ]
goodgene_exclude_target %>% dim()
# setdiff(colnames(targetgene), rownames(goodgene)) # this should be empty
othergenes <- sample_n(as.data.frame(goodgene_exclude_target), 2000) %>% 
  as.matrix() #  2000 x 1018
othergenes %>% dim

# compute pairwise distance based on correlation --------------------------
# cell.cor <- othergenes %>% cor(use="pairwise.complete.obs")
# saveRDS(cell.cor, 'data/single_cell_data/cell_cor.rds')
cell.cor <- readRDS('data/single_cell_data/cell_cor.rds')
cell.cor %>% str()
cell.dist <- as.dist(1 - abs(cell.cor))
cell.tree <- hclust(cell.dist, method="complete")
plot(cell.tree, cex=0.2)
sub_grp <- cutree(cell.tree, h=0.15)

sub_grp %>% table()


sub_grp_subset <- sub_grp[!(sub_grp %in%  which(table(sub_grp) ==1))]
sub_grp_subset %>% table()
sub_grp_subset %>% table() %>% max()
sub_grp_remove <- sub_grp[(sub_grp %in%  which(table(sub_grp) == 1))]
sub_grp_subset %>% length()
sub_grp <- sub_grp_subset


# check the clusters
cellnames <- substr(names(sub_grp), start = 0, stop = 3)
cellnames <- sub(pattern = "_$", "", x = cellnames)
clusters <- table(cellnames, sub_grp)
clusters %>% dim()
clusters %>% sum()
clusters
sub_grp %>% head(10)

sim_data <- reorder_data(sub_grp_subset = sub_grp, targetgene = targetgene)
sim_data$df %>% dim()
# saveRDS(sim_data, 'data/single_cell_data/sim_data.rds')
# sim_data <- readRDS('~/Documents/research/dag_network/data/single_cell_data/sim_data.rds')
# estimate decorrelation matrix  ------------------------------------------------------


dir.create(path = "output/single_cell23")
setwd(dir = "output/single_cell23")

saveRDS(clusters, 'clusters.rds')



networkDAG_sol_path(
  # X = res$subsetXp,
  X = sim_data$df,
  block_size=1, 
  zeropos_list = NULL,
  block_idx = sim_data$block_idx,
  lambda_len = 10,
  lambda2 = 200,
  lambda1_max_div = 50,
  maxIter = 100
)
df <- sim_data$df

Xdecor_res <- get_Xdecor(df)



# estimate CPDAG ----------------------------------------------------------

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


# test --------------------------------------------------------------------



X = sim_data$df
block_size=1 
zeropos_list = NULL
block_idx = sim_data$block_idx
lambda_len = 10
lambda2 = 100
lambda1_max_div = 50
maxIter = 100

X = X 
homega = sqrt(omega2_hat_iid)
hbeta = bhat 
htheta = thetahat
S = S 
lam1 = lambda1
lam2 = lambda2

res = 1
for(i in 1:length(block_idx)){
  tmp = det(htheta[block_idx[[i]], block_idx[[i]]])
  print(tmp)
  res = res * tmp
}

-p*log(res)
