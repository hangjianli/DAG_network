temp_runsims <- function(repp = 1, args, estimands){
  setwd(paste0("~/Dropbox/research/code/",args$setting))
  # init params ---------------------------------------------------------------------
  lambda.len = 10
  max.iter = 100
  today_date <- as.Date(args$setting, format = "%y%m%d")
  sim = repp
  X_ <- sim_X(vers = sim, p = estimands$realp,
              args = args, omg.sq = estimands$omg.sq,
              sig = estimands$sig, b = estimands$b)
  X <- X_$X
  # GES ---------------------------------------------------------------------
  n <- dim(X)[1]
  p <- dim(X)[2]
  dimnames(estimands$b) <- list(as.character(1:p), as.character(1:p))
  adjmat_trueCPDAG <- bnstruct::dag.to.cpdag(1*(estimands$b != 0)) 
  # My methods --------------------------------------------------------------
  for (sim in 1:repp){
    setwd(paste0("~/Dropbox/research/code/", args$setting,"--", sim))
    # mysims ------------------------------------------------------------------
    X_p <- Permute_X(X, seed = sim)
    XX <- X_p$Xperm
    # PC ----------------------------------------------------------------------
    BICscores_main <- readRDS('BICscores_main.rds')
    bestk_bic_main <- which.min(BICscores_main)
    bestresult_main_bic <- readRDS(file = paste0(format(today_date, "%Y-%m-%d"), "-nsim-", sim,
                                                 "-vers-", bestk_bic_main, "-lam-", bestk_bic_main,
                                                 "-mainResult",".rds"))
    XXtheta <- chol(bestresult_main_bic$theta_est) %*% XX
    res_pc <- rcausal::pc(df = XXtheta, continuous = T, depth = -1, significance = 0.01)
    adjmat_pc_CPDAG <- get_adjmat_from_pc(res_pc, estimands$realp)
    shd_pc <- compute_SHD_detail(adjmat_pc_CPDAG, adjmat_trueCPDAG)
    saveRDS(unlist(shd_pc), "shdL_pc.rds")
    cat("finished sim: ", sim, "\n")
  }
}

calc_aver <- function(args, estimands, repp = 10){
  setting <- args$setting
  total <- vector(mode = "numeric", length = 5)
  for (i in 1:repp){
    shdL_pc <- readRDS(paste0("~/Dropbox/research/code/",setting, "--", i, "/shdL_pc.rds"))  
    total <- total + shdL_pc
  }
  total <- round(total / repp, 2)
  saveRDS(total, file = paste0("../", setting, "/pcLtotal.rds"))
  return(total)
}

# 
temp_runsims(repp = 1, args = args, estimands = estimands)
calc_aver(args,estimands,repp = 5)
total <- readRDS(paste0("~/Dropbox/research/code/",args$se,"/total.rds"))
pcLtotal <- readRDS(paste0("~/Dropbox/research/code/",args$setting,"/pcLtotal.rds"))
total2 <- cbind(total, pcLtotal)
saveRDS(total2, paste0("~/Dropbox/research/code/", args$setting, "/total2.rds"))
