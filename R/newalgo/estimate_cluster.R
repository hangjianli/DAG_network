setwd("~/Documents/research/dag_network/")

targetgene <- readRDS("data/single_cell_data/sig_genes_log_val.rds")
targetgene <- t(targetgene)
targetgene %>% dim() # 51 target genes
goodgene <- readRDS('data/single_cell_data/allgoodgenes.rds')
goodgene %>% dim() # all the genes that passed VOC test
goodgene_exclude_target <- goodgene[!(rownames(goodgene) %in% colnames(targetgene)), ]
goodgene_exclude_target %>% dim()
# setdiff(colnames(targetgene), rownames(goodgene)) # this should be empty
sc_idx_full <- readRDS("~/Documents/research/dag_network/data/single_cell_data/single_cell_block_idx_full.rds")

# two-step clustering -----------------------------------------------------
# sim_data <- two_step_cluster(
#   othergenes = othergenes,
#   targetgene = targetgene,
#   sc_idx_full = sc_idx_full,
#   corr_thr = c(0.18, 0.16, 0.23, 0.2, 0.17, 0.185, 0.21)
#   )


# simple clustering -----------------------------------------------------


# compute pairwise distance based on correlation --------------------------
# cell.cor <- othergenes %>% cor(use="pairwise.complete.obs")
# saveRDS(cell.cor, 'data/single_cell_data/cell_cor.rds')
# cell.cor <- readRDS('data/single_cell_data/cell_cor.rds')
cell.cor %>% str()
cell.dist <- as.dist(1 - abs(cell.cor))
cell.tree <- hclust(cell.dist, method="complete")
plot(cell.tree, cex=0.2)
sub_grp <- cutree(cell.tree, h=0.35)

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

# sim_data <- reorder_data(sub_grp_subset = sub_grp, targetgene = targetgene)
# sim_data$df %>% dim()
# saveRDS(sim_data, 'data/single_cell_data/sim_data.rds')
# sim_data <- readRDS('~/Documents/research/dag_network/data/single_cell_data/sim_data.rds')
# estimate decorrelation matrix  ------------------------------------------------------

sim_data <- list(
  df = resdf,
  block_idx = reslist
)

# main ---------------------------------------------------------------------

sim_name = 39
dir.create(path = paste0("~/Documents/research/dag_network/output/single_cell",sim_name))
setwd(dir = paste0("~/Documents/research/dag_network/output/single_cell", sim_name))

# saveRDS(clusters, 'clusters.rds')
set.seed(1)
othergenes <- sample_n(as.data.frame(goodgene_exclude_target), 8000) %>% 
  as.matrix() #  2000 x 1018
othergenes %>% dim()
sim_data <- two_step_cluster(
  othergenes = othergenes,
  targetgene = targetgene,
  sc_idx_full = sc_idx_full,
  corr_thr = rep(0.85, 7)
)
getwd()
# sim_data <- readRDS('sim_data.rds')
saveRDS(sim_data, 'sim_data.rds')

X = apply(sim_data$df, 2, scale, scale=F)
rownames(X) <- rownames(sim_data$df)

con <- file("test.log")
sink(con, append=TRUE)
sink(con, append=TRUE, type="message")

networkDAG_sol_path(
  # X = res$subsetXp,
  X = X,
  block_size=NULL, 
  zeropos_list = NULL,
  block_idx = sim_data$block_idx,
  lambda_len = 10,
  lambda2 = 0.01,
  lambda1_max_div = 3000,
  maxIter = 100
)
sink() 
sink(type="message")

Xdecor_res <- get_Xdecor(X)

BIC_main_result <- readRDS(paste0("~/Documents/research/dag_network/output/single_cell", sim_name, "/BIC_main_result.rds"))
BIC_main_result %>% unlist()
BIC_1iter_result <- readRDS(paste0("~/Documents/research/dag_network/output/single_cell", sim_name, "/BIC_1iter_result.rds"))
BIC_1iter_result %>% unlist()
BIC_baseline_result <- readRDS(paste0("~/Documents/research/dag_network/output/single_cell", sim_name, "/BIC_baseline_result.rds"))
BIC_baseline_result %>% unlist()
# estimate CPDAG ----------------------------------------------------------

GES_sol(X, block_idx = sim_data$block_idx, decor = F)
GES_sol(Xdecor_res$X_decor, block_idx = sim_data$block_idx, decor = T)
# GES_sol(Xdecor_res$X_decor_1iter, decor = T)
pc_sol(X, block_idx = sim_data$block_idx, decor = F)
pc_sol(Xdecor_res$X_decor, block_idx = sim_data$block_idx, decor = T)
sparsebn_sol(X, block_idx = sim_data$block_idx, decor = F)
sparsebn_sol(Xdecor_res$X_decor, block_idx = sim_data$block_idx, decor = T)

bic_baseline <- readRDS("BIC_baseline_result.rds")
bic_1iter <- readRDS("BIC_1iter_result.rds")
bic_main <- readRDS("BIC_main_result.rds")
bic_ges <- readRDS("fGES_BIC_result.rds")
bic_pc <- readRDS("pc_BIC_result.rds")
bic_sbn <- readRDS("sbn_BIC_result.rds")
bic_ges_decor <- readRDS("fGES_BIC_result_decor.rds")
bic_pc_decor <- readRDS("pc_BIC_result_decor.rds")
bic_sbn_decor <- readRDS("sbn_BIC_result_decor.rds")

all_bic <- data.frame(
  bic_baseline = bic_baseline %>% unlist(),
  bic_1iter = bic_1iter %>% unlist(),
  bic_main = bic_main %>% unlist(),
  bic_ges = bic_ges %>% unlist(),
  bic_ges_decor = bic_ges_decor %>% unlist(),
  bic_pc = bic_pc %>% unlist(),
  bic_pc_decor = bic_pc_decor %>% unlist(),
  bic_sbn = bic_sbn %>% unlist(),
  bic_sbn_decor = bic_sbn_decor %>% unlist()
)
saveRDS(all_bic, file = "all_bic.rds")
all_bic


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



X = apply(sim_data$df, 2, scale, scale=F)
rownames(X) <- rownames(sim_data$df)

block_size=NULL
zeropos_list = NULL
block_idx = sim_data$block_idx
lambda_len = 10
lambda2 = 200
lambda1_max_div = 3000
maxIter = 100

k=8
lambda1 = lambda.path[k]
tol = 1e-7


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

