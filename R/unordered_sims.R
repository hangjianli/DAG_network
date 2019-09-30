unordered_runsims <- function(start_pos = 1, repp = 1, args, estimands, thr, 
                              max.iter = 100, lambda.len = 10, div = 100){
  # dir.create(path = args$setting)
  setwd(paste0("~/Dropbox/research/code/",args$setting))
  today_date <- Sys.Date()
  
  # My methods --------------------------------------------------------------
  for (sim in start_pos:repp){
    dir.create(path = paste0("~/Dropbox/research/code/", args$setting,"--", sim))
    setwd(paste0("~/Dropbox/research/code/", args$setting,"--", sim))
    
    estimands <- generate_parameters(args, seed = sim*2, theta_name = estimands$theta_name, 
                                     bname = estimands$bname, btype = estimands$btype)
    saveRDS(estimands, file = "estimands.rds")
    # sim new X ---------------------------------------------------------------
    X_ <- sim_X(vers = sim, p = estimands$realp,
                args = args, omg.sq = estimands$omg.sq,
                sig = estimands$sig, b = estimands$b)
    X <- X_$X
    saveRDS(X_, file = paste0(format(today_date, "%Y-%m-%d"), "-vers-",
                              args$setting, "-X-",sim,".rds"))
    # GES ---------------------------------------------------------------------
    n <- dim(X)[1]
    p <- dim(X)[2]
    dimnames(estimands$b) <- list(as.character(1:p), as.character(1:p))
    adjmat_trueCPDAG <- bnstruct::dag.to.cpdag(1*(estimands$b != 0)) 
    saveRDS(adjmat_trueCPDAG, "adjmat_trueCPDAG.rds")
    set.seed(100)
    fgsX <- fges(df = X, penaltydiscount = 2.0, maxDegree = -1,
                 faithfulnessAssumed = TRUE, verbose = F)
    adjmat_fgesCPDAG_X <- get_adjmat_from_fges(fgsX$edges, p = p, varnames = dimnames(X)[[2]])
    shdX_ges <- unlist(compute_SHD_detail(adjmat_fgesCPDAG_X, adjmat_trueCPDAG, estimands$s0))
    saveRDS(shdX_ges, paste0("X-",sim,"-shdX_ges.rds"))
    
    
    # mysims ------------------------------------------------------------------
    X_p <- Permute_X(X, seed = sim)
    saveRDS(X_p, file = paste0(format(today_date, "%Y-%m-%d"), "-vers-", sim, "-permutation",".rds"))
    XX <- X_p$Xperm
    rho.est <- rep(1,p)
    # PC ----------------------------------------------------------------------
    res_pc <- rcausal::pc(df = XX, continuous = T, depth = -1, significance = 0.001)
    # res_pc <- pcalg::pc(suffStat = XX, indepTest = gaussCItest, alpha = 0.001)
    adjmat_pc_CPDAG <- get_adjmat_from_pc(res_pc, estimands$realp)
    shd_pc <- unlist(compute_SHD_detail(adjmat_pc_CPDAG, adjmat_trueCPDAG, estimands$s0))
    saveRDS(shd_pc, "shd_pc.rds")
    # get lambda path ---------------------------------------------------------
    lambda.path <- get_lam_path(p, XX, rho.est, lambda.len, div = div)
    saveRDS(lambda.path, file = "lambda_path.rds")
    # some initializaiton -----------------------------------------------------
    # BICscores_main <- BICscores_flipflop <- BICscores_bench <- rep(0, lambda.len)
    # minrowcor_main <- minrowcor_bench <- minrowcor_flip <- rep(0,lambda.len)
    BICscores_main <- minrowcor_main <- rep(0, lambda.len)
    s0_seq <- rep(0, lambda.len)
    for(k in 1:lambda.len){
      set.seed(1)
      lambda <- rev(lambda.path)[k]
      # Main algorithm ----------------------------------------------------------
      result_main <- main_iteration_fast2(X = XX, vers = k*2, lambda = lambda,
                                         lambda2 = 0.01*n, max.iter = max.iter,
                                         estimands = estimands, args = args)
      
      s0_seq[k] <- result_main$s0_est
      if(any(apply(result_main$b_est, 2, function(x) sum(abs(x) > 1e-4)) >= n)){ #prevent high dimensional MLE
        cat("[INFO] lambda is too small, algorithm stops at k = ", k, "\n")
        break
      }
      # save main result
      saveRDS(result_main, file = paste0(format(today_date, "%Y-%m-%d"), "-nsim-", sim,"-vers-",
                                         k, "-lam-", k, "-mainResult",".rds"))
      # save plots
      png(paste0("NLikelihood", k, "-sim", ".png"))
      plot(result_main$lik_seq, pch = 16, type = "b", cex = 0.7, ylab = "Negative LL", xlab = "Iter")
      dev.off()
      # # --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
      cor_est <- cor(t(chol(result_main$theta_est)%*%(XX - XX%*%result_main$b_est)))
      minrowcor_main[k] <- sum(abs(cor_est[upper.tri(cor_est)]))
      mle_result <- dag_mle_estimation(X = XX, Bhat = result_main$b_est, 
                                       Lhat = chol(result_main$theta_est))
      # save MLE result from main
      saveRDS(mle_result, file = paste0(format(today_date, "%Y-%m-%d"),  "-nsim-", sim,
                                        "-vers-", k, "-lam-", k,"-mainMLEResult",".rds"))
      BIC_result <- BIC_dag(X = XX, bmle = mle_result$Bmle, omgmle = mle_result$omgmlesq,
                            theta = result_main$theta_est)
      BICscores_main[k] <- BIC_result$BIC
      saveRDS(BICscores_main, "BICscores_main.rds")
      saveRDS(minrowcor_main, "minrowcor_main.rds")
      # benchmark ---------------------------------------------------------------
      # result_bench <- lassoIdentTheta(X = XX, lambda = lambda, vers = k*2, zeropos = estimands$zeropos,
      #                                 max.iter = 50, lambda2 = 0.01)
      # saveRDS(result_bench, file = paste0(format(today_date, "%Y-%m-%d"), "-nsim-", sim,"-vers-",
      #                                     k, "-lam-", k, "-benchResult",".rds"))
      # save plots
      # png(paste0("NLikelihood_bench", k, "-sim", ".png"))
      # plot(result_bench$likeli_seq, pch = 16, type = "b", cex = 0.7, ylab = "Negative LL", xlab = "Iter")
      # dev.off()
      # --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
      # cor_est <- cor(XX - XX%*%result_bench$b_est)
      # minrowcor_bench[k] <- sum(abs(cor_est[upper.tri(cor_est)]))
      # mle_result_bench <- dag_mle_estimation(X = XX, Bhat = result_bench$b_est, Lhat = diag(n))
      # saveRDS(mle_result_bench, file = paste0(format(today_date, "%Y-%m-%d"),  "-nsim-", sim,
      #                                         "-vers-", k, "-lam-", k,"-benchMLEResult",".rds"))
      # BIC_result <- BIC_dag(X = XX, bmle = mle_result_bench$Bmle,
      #                       omgmle = mle_result_bench$omgmlesq, theta = diag(n))
      # BICscores_bench[k] <- BIC_result$BIC
      # saveRDS(BICscores_bench, "BICscores_bench.rds")
      # saveRDS(minrowcor_bench, "minrowcor_bench.rds")

      cat("[INFO] lambda ", k, "is finished.\n")
    }

    # saveRDS(BICscores_bench[BICscores_bench!=0], "BICscores_bench.rds")
    saveRDS(BICscores_main[BICscores_main!=0], "BICscores_main.rds")
    # saveRDS(minrowcor_bench[minrowcor_bench!=0], "minrowcor_bench.rds")
    saveRDS(minrowcor_main[minrowcor_main!=0],"minrowcor_main.rds")
    
    
    png(paste0("BICscores_main", sim ,".png"))
    plot(BICscores_main[BICscores_main!=0], pch = 16, type = "b")
    dev.off()
    png(paste0("minrowcor_main", sim ,".png"))
    plot(minrowcor_main[minrowcor_main!=0], pch = 16, type = "b")
    dev.off()
    # png(paste0("BICscores_bench", sim ,".png"))
    # plot(BICscores_bench, pch = 16, type = "b")
    # dev.off()
    # png(paste0("minrowcor_bench", sim ,".png"))
    # plot(minrowcor_bench, pch = 16, type = "b")
    # dev.off()
    
    
    bestk_bic_main <- which.min(BICscores_main[BICscores_main!=0])
    bestk_cor_main <- which.min(minrowcor_main[minrowcor_main!=0])
    # bestk_bic_bench <- which.min(BICscores_bench)
    # bestk_cor_bench <- which.min(minrowcor_bench)
    
    
    allShdS <- get_shd_for_several_methods(kmainbic = bestk_bic_main,
                                           # kbenchbic = bestk_bic_bench,
                                           kmaincor = bestk_cor_main,
                                           # kbenchcor = bestk_cor_bench,
                                           today_date = today_date,
                                           tt = sim,
                                           gesshd = shdX_ges, 
                                           pcshd = shd_pc,
                                           adjmat_trueCPDAG = adjmat_trueCPDAG,
                                           thresh = thr,
                                           estimands = estimands,
                                           ordered = F)
    
    
    
    bestresult_main_bic <- readRDS(file = paste0(format(today_date, "%Y-%m-%d"), "-nsim-", sim,
                                                 "-vers-", bestk_bic_main, "-lam-", bestk_bic_main,
                                                 "-mainResult",".rds"))
    thetaRE <- norm(bestresult_main_bic$theta_est - estimands$theta, "f") / norm(estimands$theta, "f")
    # image(as(round(bestresult_main_bic$theta_est,3), class(estimands$theta)))
    # image(estimands$theta)
    saveRDS(thetaRE, file = "thetaRE.rds")
    

    # GESL --------------------------------------------------------------------
    XXtheta <- chol(bestresult_main_bic$theta_est) %*% XX
    dimnames(XXtheta) <- dimnames(XX)
    fgs_theta <- fges(df = XXtheta, penaltydiscount = 2.0, maxDegree = -1,
                      faithfulnessAssumed = TRUE, verbose = F)
    adjmat_fgesCPDAG_theta <- get_adjmat_from_fges(fgs_theta$edges, p = p,
                                                   varnames = dimnames(X)[[2]])
    shdXLBIC <- compute_SHD_detail(adjmat_fgesCPDAG_theta, adjmat_trueCPDAG, estimands$s0)
    saveRDS(unlist(shdXLBIC), "shdXLBIC.rds")

    # PCL ---------------------------------------------------------------------
    res_pc <- rcausal::pc(df = XXtheta, continuous = T, depth = -1, significance = 0.01)
    adjmat_pc_CPDAG <- get_adjmat_from_pc(res_pc, p = p)
    shdL_pc <- unlist(compute_SHD_detail(adjmat_pc_CPDAG, adjmat_trueCPDAG, estimands$s0))
    saveRDS(shdL_pc, "shdL_pc.rds")
    
    
    allShdS$shdXLBIC <- unlist(shdXLBIC)
    allShdS$shdL_pc <- unlist(shdL_pc) 
    
    saveRDS(data.frame(allShdS), file = "allshd.rds")
    if(bestk_bic_main != bestk_cor_main){
      bestresult_main_cor <- readRDS(file = paste0(format(today_date, "%Y-%m-%d"), "-nsim-", sim,
                                                   "-vers-", bestk_cor_main, "-lam-", bestk_cor_main,
                                                   "-mainResult",".rds"))
      XXtheta <- chol(bestresult_main_cor$theta_est) %*% XX
      dimnames(XXtheta) <- dimnames(XX)
      fgs_theta <- fges(df = XXtheta, penaltydiscount = 2.0, maxDegree = -1,
                        faithfulnessAssumed = TRUE, verbose = F)
      adjmat_fgesCPDAG_theta <- get_adjmat_from_fges(fgs_theta$edges, p = p,
                                                     varnames = dimnames(X)[[2]])
      shdXXLcor <- compute_SHD_detail(adjmat_fgesCPDAG_theta, adjmat_trueCPDAG, estimands$s0)
      saveRDS(shdXXLcor, "shdXXLcor.rds")
    }
    cat("---------------------------------------------------------------\n")
    cat("[INFO] lambda ", k, "is finished.\n")
    cat("---------------------------------------------------------------\n")
  }
}
