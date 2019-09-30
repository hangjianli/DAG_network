ordered_runsims <- function(start = 1, repp = 5, args){
  
  setwd(paste0("~/Dropbox/research/code/",args$setting))
  lambda.len = 10
  max.iter = 100
  today_date <- Sys.Date()
  
  for (sim in start:repp){
    dir.create(path = paste0("~/Dropbox/research/code/",args$setting,"--", sim))
    setwd(paste0("~/Dropbox/research/code/",args$setting,"--", sim))
    # generate a new dag ------------------------------------------------------
    estimands <- generate_parameters(args, seed = sim, theta_name = estimands$theta_name, 
                                     bname = estimands$bname, btype = estimands$btype)
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
    
    zeropos_psi= NULL
    dimnames(estimands$b) <- list(as.character(1:p), as.character(1:p))
    adjmat_trueCPDAG <- bnstruct::dag.to.cpdag(1*(estimands$b != 0)) 
    XX <- X
    # get lambda path ---------------------------------------------------------
    rho.est <- rep(1,p)
    lambda.path <- get_lam_path(p, XX, rho.est, lambda.len, 1000)
    saveRDS(lambda.path, file = "lambda_path.rds")
    # some initializaiton -----------------------------------------------------
    BICscores_main <- BICscores_flipflop <- BICscores_bench <- rep(0, lambda.len)
    minrowcor_main <- minrowcor_bench <- minrowcor_flip <- rep(0,lambda.len)
    thr <- 0
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
      
      png(paste0("NLikelihood", k, "-sim", ".png"))
      plot(result_main$lik_seq, pch = 16, type = "b", cex = 0.7, ylab = "Negative LL", xlab = "Iter")
      dev.off()
      # --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
      cor_est <- cor(t(chol(result_main$theta_est)%*%(XX - XX%*%result_main$b_est)))
      minrowcor_main[k] <- sum(abs(cor_est[upper.tri(cor_est)]))
      mle_result <- dag_mle_estimation(X = XX, Bhat = result_main$b_est,
                                       Lhat = chol(result_main$theta_est))
      thr <- max(thr, mle_result$threshold)
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
      # --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
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

      # flipflop ----------------------------------------------------------------
      set.seed(1)
      resultflip <- flipflop(X = XX, n = args$n,
                             args = args,
                             zeropos_psi = ,
                             zeropos_list = estimands$zeropos_list,
                             p = estimands$realp,
                             max.iter = 10,
                             lambda = lambda,
                             lambda2 = 0.1*n)
      #
      psi_hat <- resultflip$psi_est
      theta_hat <- resultflip$theta_est
      psi_zeros <- which(abs(psi_hat) < 1e-4, arr.ind = T)
      # theta_zeros <- which(abs(theta_hat) < 1e-4, arr.ind = T)
      theta_zeros <- get_zeros(theta_hat,args,estimands)
      saveRDS(resultflip, file = paste0(format(today_date, "%Y-%m-%d"),
                                       "-nsim-", sim,
                                       "-vers-", k, "-lam-",
                                       k,"-flipflop",".rds"))

      # if(any(diff(resultflip$likeli_seq[resultflip$likeli_seq!=0]) > 0))
      #   warning("likelihood not descreasing ! \n")

      png(paste0("B_est_flipflop_lik_", k ,".png"))
      plot(resultflip$likeli_seq[resultflip$likeli_seq!=0], pch = 16, type = "b")
      lines(resultflip$likeli_seq[resultflip$likeli_seq!=0], type = "o")
      dev.off()

      cat("[INFO] fliplop MLE calculation. \n")
      resultMLE <- flipflop(X=XX, n = n, p = p,
                            args = args, zeropos_list = theta_zeros,
                            zeropos_psi = psi_zeros,
                            max.iter = 10,
                            lambda = 0.001,
                            lambda2 = 0.001)
      saveRDS(resultMLE, file = paste0(format(today_date, "%Y-%m-%d"),
                                       "-nsim-", sim,
                                       "-vers-", k, "-lam-",
                                       k,"-flipflopMLE",".rds"))
      
      # resultMLE <- readRDS(file = paste0(format(today_date, "%Y-%m-%d"), 
      #                                    "-nsim-", sim,
      #                                    "-vers-", k, "-lam-", 
      #                                    k,"-flipflopMLE",".rds"))
      
      BICscores_flipflop[k] <- -n*log(det(resultMLE$psi_est)) - p*log(det(resultMLE$theta_est)) +
        sum(diag(resultMLE$theta_est%*%XX%*%resultMLE$psi_est%*%t(XX))) + log(max(n,p))*sum(abs(resultMLE$psi_est) > 1e-4)
      cor_est <- cor(t(chol(resultMLE$theta_est)%*%(XX - XX%*%resultMLE$b_est)))
      minrowcor_flip[k] <- sum(abs(cor_est[upper.tri(cor_est)]))
      saveRDS(minrowcor_flip,"minrowcor_flip.rds")
      saveRDS(BICscores_flipflop,"BICscores_flipflop.rds")
      cat("[INFO] lambda ", k, "is finished.\n")
    }
    saveRDS(BICscores_bench[BICscores_bench!=0], "BICscores_bench.rds")
    saveRDS(BICscores_main[BICscores_main!=0], "BICscores_main.rds")
    saveRDS(minrowcor_bench[minrowcor_bench!=0], "minrowcor_bench.rds")
    saveRDS(minrowcor_main[minrowcor_main!=0],"minrowcor_main.rds")
    saveRDS(BICscores_flipflop[BICscores_flipflop!=0],"BICscores_flipflop.rds")
    saveRDS(minrowcor_flip[minrowcor_flip!=0], "minrowcor_flip.rds")
    
    
    png(paste0("BICscores_flipflop", sim ,".png"))
    plot(BICscores_flipflop[BICscores_flipflop!=0], pch = 16, type = "b")
    dev.off()
    png(paste0("minrowcor_flip", sim ,".png"))
    plot(minrowcor_flip[minrowcor_flip!=0], pch = 16, type = "b")
    dev.off()
    
    
    png(paste0("BICscores_main", sim ,".png"))
    plot(BICscores_main[BICscores_main!=0], pch = 16, type = "b")
    dev.off()
    png(paste0("minrowcor_main", sim ,".png"))
    plot(minrowcor_main[minrowcor_main!=0], pch = 16, type = "b")
    dev.off()
    png(paste0("BICscores_bench", sim ,".png"))
    plot(BICscores_bench, pch = 16, type = "b")
    dev.off()
    png(paste0("minrowcor_bench", sim ,".png"))
    plot(minrowcor_bench, pch = 16, type = "b")
    dev.off()
    bestk_bic_flop <- which.min(BICscores_flipflop)
    bestk_cor_flop <- which.min(minrowcor_flip)
   
    bestk_bic_main <- which.min(BICscores_main)
    bestk_cor_main <- which.min(minrowcor_main)
    bestk_bic_bench <- which.min(BICscores_bench)
    bestk_cor_bench <- which.min(minrowcor_bench)
    bestk_bic_flip <- which.min(BICscores_flipflop)
    bestk_cor_flip <- which.min(minrowcor_flip)
    
    allShdS <- get_shd_for_several_methods(kmainbic = bestk_bic_main,
                                           kbenchbic = bestk_bic_bench,
                                           kmaincor = bestk_cor_main,
                                           kbenchcor = bestk_cor_bench,
                                           kflipbic = bestk_bic_flip,
                                           kflipcor = bestk_cor_flip,
                                           today_date = today_date,
                                           tt = sim,
                                           gesshd = NULL, 
                                           pcshd = NULL,
                                           adjmat_trueCPDAG = adjmat_trueCPDAG,
                                           thresh = 0.2,
                                           estimands = estimands)
    
    saveRDS(allShdS, file = "allshd.rds")
    # bestresult_main_bic <- readRDS(file = paste0(format(today_date, "%Y-%m-%d"), "-nsim-", sim,
                                                 # "-vers-", bestk_bic_main, "-lam-", bestk_bic_main,
                                                 # "-mainResult",".rds"))
    # thetaRE <- norm(bestresult_main_bic$theta_est - estimands$theta, "f") / norm(estimands$theta, "f")
    # saveRDS(thetaRE, file = "thetaRE.rds")
    
    
  }
}
