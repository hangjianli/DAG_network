

solution2_idx <- select.parameter(test2, XX)
solution2 <- select(test2, index = solution2_idx)
myadj <- ifelse(abs(estimands$b) > 0, 1, 0)
sparsebngraph <- graph_from_adjacency_matrix(myadj, mode = "directed")

plot(sparsebngraph, vertex.label = NA,
     vertex.size = 5,
     vertex.label.color = gray(0),
     vertex.color = gray(0.9),
     edge.color = gray(0),
     edge.arrow.size = 0.15)

plotDAG(solution2)
num.edges(solution2)

newX <- sparsebnData(XXtheta, type = "continuous")
test2 <- estimate.dag(newX)
solution2_idx <- select.parameter(test2, newX)
solution2 <- select(test2, index = solution2_idx)

sparsebn_origin <- get.adjacency.matrix(solution)

compute_SHD_detail(sparsebn_origin, adjmat_trueCPDAG, estimands$s0)



add_sparsebn <- function(start = 1, sim_num = 10, setting, today_date, seed=100){
  for(sim in start:sim_num){
    cat(paste0("[INFO] start sim ", sim, " !\n"))
    setwd(paste0("~/Dropbox/research/code/", setting,"--", sim))
    estimands <- readRDS("estimands.rds")
    X_ <- readRDS(paste0(format(today_date, "%Y-%m-%d"), "-vers-",
                         setting, "-X-",sim,".rds"))  
    X <- X_$X
    
    n <- dim(X)[1]
    p <- dim(X)[2]
    dimnames(estimands$b) <- list(as.character(1:p), as.character(1:p))
    adjmat_trueCPDAG <- readRDS("adjmat_trueCPDAG.rds")
    
    X_p <- readRDS(paste0(format(today_date, "%Y-%m-%d"), "-vers-", sim, "-permutation",".rds"))
    XX <- X_p$Xperm
    
    set.seed(seed)
    sbX <- sparsebnData(XX, type = "continuous")
    sb_path <- estimate.dag(sbX)
    sol_idx <- select.parameter(sb_path, sbX)
    sol_base <- select(sb_path, index = sol_idx)
    adjmat_sbbase <- get.adjacency.matrix(sol_base)
    shdX_sb <- unlist(compute_SHD_detail(adjmat_sbbase, adjmat_trueCPDAG, estimands$s0))
    saveRDS(shdX_sb, paste0("X-",sim,"-shdX_sb.rds"))
    
    BICscores_main <- readRDS("BICscores_main.rds")
    minrowcor_main <- readRDS("minrowcor_main.rds")
    bestk_bic_main <- which.min(BICscores_main[BICscores_main!=0])
    bestk_cor_main <- which.min(minrowcor_main[minrowcor_main!=0])
    bestresult_main_bic <- readRDS(file = paste0(format(today_date, "%Y-%m-%d"), "-nsim-", sim,
                                                 "-vers-", bestk_bic_main, "-lam-", bestk_bic_main,
                                                 "-mainResult",".rds"))
    
    XXtheta <- chol(bestresult_main_bic$theta_est) %*% XX
    dimnames(XXtheta) <- dimnames(XX)
    sbXtheta <- sparsebnData(XXtheta, type = "continuous")
    sb_path_theta <- estimate.dag(sbXtheta)
    sol_idx_theta <- select.parameter(sb_path_theta, sbXtheta)
    sol_base_deco <- select(sb_path_theta, index = sol_idx_theta)
    adjmat_sbbase_decor <- get.adjacency.matrix(sol_base_deco)
    shdX_sb_decor <- unlist(compute_SHD_detail(adjmat_sbbase_decor, adjmat_trueCPDAG, estimands$s0))
    saveRDS(shdX_sb_decor, paste0("X-",sim,"-shdX_sb_decor.rds"))
    allshd <- readRDS("allshd.rds")
    allshd$shdSB <- shdX_sb
    allshd$shdSBL <- shdX_sb_decor
    saveRDS(allshd, "allshd.rds")
    
    cat(paste0("[INFO] sim ", sim, " done.\n"))
  }
}


