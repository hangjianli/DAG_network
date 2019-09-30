ordered_runsims_new <- function(start = 1, repp = 5, args, estimands,
                                max.iter = 100, lam_div=100, lambda.len = 10){
  ####################################
  # 4/9/2019
  # In this version, when n > p, flipflop is not performed on lambda1 path but given the true support.
  #
  ####################################
  setwd(paste0("~/Dropbox/research/code/",args$setting))
  today_date <- Sys.Date()
  
  for (sim in start:repp){
    dir.create(path = paste0("~/Dropbox/research/code/",args$setting,"--", sim))
    setwd(paste0("~/Dropbox/research/code/",args$setting,"--", sim))
    # generate a new dag ------------------------------------------------------
    estimands <- generate_parameters(args, seed = sim*3, theta_name = estimands$theta_name, 
                                     bname = estimands$bname, btype = estimands$btype)
    estimands$sig <- diag(args$n)
    estimands$theta <- diag(args$n)
    saveRDS(estimands, file = "estimands.rds")
    # generate data X ---------------------------------------------------------
    X_ <- sim_X(vers = sim, p = estimands$realp,
                args = args, omg.sq = estimands$omg.sq,
                sig = estimands$sig, b = estimands$b)
    X <- X_$X
    saveRDS(X_, file = paste0(format(today_date, "%Y-%m-%d"), "-vers-",
                              args$setting, "-X-",sim,".rds"))
    n <- dim(X)[1]
    p <- dim(X)[2]
    
    zeropos_psi_inv= NULL
    dimnames(estimands$b) <- list(as.character(1:p), as.character(1:p))
    # adjmat_trueCPDAG <- bnstruct::dag.to.cpdag(1*(estimands$b != 0)) 
    adjmat_trueCPDAG <- 1*(estimands$b != 0)
    dimnames(adjmat_trueCPDAG) <- dimnames(estimands$b)
    saveRDS(adjmat_trueCPDAG, "adjmat_true_DAG.rds")
    XX <- X
    # compute_psi and flipflop-------------------------------------------------------------
    psi_true_inv <- t(diag(p) - estimands$b) %*% diag(1/estimands$omg.sq)%*% t(diag(p) - estimands$b)
    zeropos_psi_inv <- which(abs(psi_true_inv) < 1e-8, arr.ind = T)
    # dim(zeropos_psi_inv)
    # set.seed(1)
    # resultflip <- flipflop(X = XX, n = n,
    #                        args = args,
    #                        zeropos_psi = zeropos_psi_inv,
    #                        zeropos_list = estimands$zeropos_list,
    #                        p = estimands$realp,
    #                        max.iter = 20,
    #                        lambda = 0.0001,
    #                        lambda2 = 0.0001)
    # saveRDS(resultflip, file = paste0(format(today_date, "%Y-%m-%d"),
    #                                   "-nsim-", sim,
    #                                   "-flipflop",".rds"))
    # some initializaiton -----------------------------------------------------
    rho.est <- rep(1,p)
    lambda.path <- get_lam_path(p, XX, rho.est, lambda.len, lam_div)
    saveRDS(lambda.path, file = "lambda_path.rds")
    BICscores_main <- BICscores_flipflop <- BICscores_bench <- rep(0, lambda.len)
    minrowcor_main <- minrowcor_bench <- minrowcor_flip <- rep(0,lambda.len)

    # Iterative through lambda ------------------------------------------------
    for(k in 1:lambda.len){
      set.seed(1)
      lambda <- rev(lambda.path)[k]
      # Main algorithm ----------------------------------------------------------
      result_main <- main_iteration_fast2(X = XX, vers = k*2, lambda = lambda,
                                          lambda2 = 0.01, max.iter = max.iter,
                                          estimands = estimands, args = args)
      
      if(any(apply(result_main$b_est, 2, function(x) sum(abs(x) > 1e-4)) >= n)){
        cat("[INFO] lambda is too small, algorithm stops at k = ", k, "\n")
        break
      }
      saveRDS(result_main, file = paste0(format(today_date, "%Y-%m-%d"), "-nsim-", sim,"-vers-",
                                         k, "-lam-", k, "-mainResult",".rds"))
      
      # result_main <- readRDS(file = paste0(format(today_date, "%Y-%m-%d"), "-nsim-", sim,"-vers-",
      #                                             k, "-lam-", k, "-mainResult",".rds"))
      
      # png(paste0("NLikelihood", k, "-sim", ".png"))
      # plot(result_main$lik_seq, pch = 16, type = "b", cex = 0.7, ylab = "Negative LL", xlab = "Iter")
      # dev.off()
      # --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
      cor_est <- cor(t(chol(result_main$theta_est)%*%(XX - XX%*%result_main$b_est)))
      minrowcor_main[k] <- sum(abs(cor_est[upper.tri(cor_est)]))
      mle_result <- dag_mle_estimation(X = XX, Bhat = result_main$b_est,
                                       Lhat = chol(result_main$theta_est))
      saveRDS(mle_result, file = paste0(format(today_date, "%Y-%m-%d"),  "-nsim-", sim,
                                        "-vers-", k, "-lam-", k,"-mainMLEResult",".rds"))
      # mle_result <- readRDS(file = paste0(format(today_date, "%Y-%m-%d"),  "-nsim-", sim,
      # "-vers-", k, "-lam-", k,"-mainMLEResult",".rds"))
      
      BIC_result <- BIC_dag(X = XX, bmle = mle_result$Bmle, omgmle = mle_result$omgmlesq,
                            theta = result_main$theta_est)
      BICscores_main[k] <- BIC_result$BIC
      saveRDS(BICscores_main, "BICscores_main.rds")
      saveRDS(minrowcor_main, "minrowcor_main.rds")
      # benchmark ---------------------------------------------------------------
      result_bench <- lassoIdentTheta(X = XX, lambda = lambda, vers = k*2,
                                      zeropos = estimands$zeropos,
                                      max.iter = 50, lambda2 = 0.01)
      
      saveRDS(result_bench, file = paste0(format(today_date, "%Y-%m-%d"),
                                          "-nsim-", sim,"-vers-",
                                          k, "-lam-", k, "-benchResult",".rds"))
      # save plots
      png(paste0("NLikelihood_bench", k, "-sim", ".png"))
      plot(result_bench$likeli_seq, pch = 16, type = "b", cex = 0.7, ylab = "Negative LL", xlab = "Iter")
      dev.off()
      # MLE of bench --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
      # result_bench <- readRDS(file = paste0(format(today_date, "%Y-%m-%d"),
      #                                       "-nsim-", sim,"-vers-",
      #                                       k, "-lam-", k, "-benchResult",".rds"))
      cor_est <- cor(XX - XX%*%result_bench$b_est)
      minrowcor_bench[k] <- sum(abs(cor_est[upper.tri(cor_est)]))
      mle_result_bench <- dag_mle_estimation(X = XX, Bhat = result_bench$b_est, Lhat = diag(n))
      saveRDS(mle_result_bench, file = paste0(format(today_date, "%Y-%m-%d"),  "-nsim-", sim,
                                              "-vers-", k, "-lam-", k,"-benchMLEResult",".rds"))
      # mle_result_bench <- readRDS(file = paste0(format(today_date, "%Y-%m-%d"),  "-nsim-", sim,
      # "-vers-", k, "-lam-", k,"-benchMLEResult",".rds"))
      BIC_result <- BIC_dag(X = XX, bmle = mle_result_bench$Bmle, 
                            omgmle = mle_result_bench$omgmlesq, theta = diag(n))
      BICscores_bench[k] <- BIC_result$BIC
      saveRDS(BICscores_bench, "BICscores_bench.rds")
      saveRDS(minrowcor_bench, "minrowcor_bench.rds")
      
      cat("[INFO] lambda ", k, "is finished.\n")
    }
    
    saveRDS(BICscores_bench[BICscores_bench!=0], "BICscores_bench.rds")
    saveRDS(BICscores_main[BICscores_main!=0], "BICscores_main.rds")
    saveRDS(minrowcor_bench[minrowcor_bench!=0], "minrowcor_bench.rds")
    saveRDS(minrowcor_main[minrowcor_main!=0],"minrowcor_main.rds")
    
    
    png(paste0("BICscores_main", sim ,".png"))
    plot(BICscores_main[BICscores_main!=0], pch = 16, type = "b")
    dev.off()
    png(paste0("minrowcor_main", sim ,".png"))
    plot(minrowcor_main[minrowcor_main!=0], pch = 16, type = "b")
    dev.off()
    png(paste0("BICscores_bench", sim ,".png"))
    plot(BICscores_bench[BICscores_bench!=0], pch = 16, type = "b")
    dev.off()
    png(paste0("minrowcor_bench", sim ,".png"))
    plot(minrowcor_bench[minrowcor_bench!=0], pch = 16, type = "b")
    dev.off()
    
    bestk_bic_main <- which.min(BICscores_main[BICscores_main!=0])
    bestk_cor_main <- which.min(minrowcor_main[minrowcor_main!=0])
    bestk_bic_bench <- which.min(BICscores_bench[BICscores_bench!=0])
    bestk_cor_bench <- which.min(minrowcor_bench[minrowcor_bench!=0])
    
    
    allShdS <- get_shd_for_several_methods(kmainbic = bestk_bic_main,
                                           kbenchbic = bestk_bic_bench,
                                           kmaincor = bestk_cor_main,
                                           kbenchcor = bestk_cor_bench,
                                           today_date = today_date,
                                           tt = sim,
                                           gesshd = NULL, 
                                           pcshd = NULL,
                                           adjmat_trueCPDAG = adjmat_trueCPDAG,
                                           thresh = 0.1,
                                           estimands = estimands,
                                           ordered = T)
    
    saveRDS(allShdS, file = "allshd2.rds")
  }
}


