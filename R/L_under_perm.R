setwd("~/Documents/research/dag_network/R")
args <- args_for_parameter()
dir.create(path = args$setting)
setwd(args$setting)
estimands <- generate_parameters(args, seed = 2)

# Plot B ------------------------------------------------------------------
# colfunc <- colorRampPalette(c("black", "white"))
png("TrueB.png")
heatmap.2(abs(estimands$b), 
          # col = colfunc(2),
          dendrogram = "none", Rowv = F, Colv = F, 
          trace = "none")
dev.off()

png("TrueTheta.png")
heatmap.2(abs(estimands$theta), 
          dendrogram = "none", 
          Rowv = F, Colv = F, 
          trace = "none")
dev.off()

png("TrueSigma.png")
heatmap.2(abs(estimands$sig), 
          dendrogram = "none", 
          Rowv = F, Colv = F, 
          trace = "none")
dev.off()
today_date <- Sys.Date()
saveRDS(estimands, file = paste0(format(today_date, "%Y-%m-%d"), "-vers-",
                                 args$setting, "-estimands",".rds"))
saveRDS(args, file = paste0(format(today_date, "%Y-%m-%d"), "-vers-", 
                            args$setting, "-args",".rds"))

# Simulate X --------------------------------------------------------------
X_ <- sim_X(vers = 1, p = estimands$realp,
            args = args, omg.sq = estimands$omg.sq,
            sig = estimands$sig, b = estimands$b)

X <- X_$X
saveRDS(X_, file = paste0(format(today_date, "%Y-%m-%d"), "-vers-",
                          args$setting, "-X",".rds"))
# run a number of permutes ------------------------------------------------
nsim = 5
lambda.len=10
max.iter = 100
tt=1

SHD_compare <- function(X, estimands, nsim = 5, lambda.len=10, max.iter = 50){
  n <- dim(X)[1]
  p <- dim(X)[2]
  dimnames(estimands$b) <- list(as.character(1:p), as.character(1:p))
  adjmat_trueCPDAG <- bnstruct::dag.to.cpdag(1*(estimands$b != 0))
  set.seed(100)
  fgsX <- fges(df = X, penaltydiscount = 2.0, maxDegree = -1,
               faithfulnessAssumed = TRUE, verbose = F)
  adjmat_fgesCPDAG_X <- get_adjmat_from_fges(fgsX$edges, p = p, varnames = dimnames(X)[[2]])
  shdX <- unlist(compute_SHD_detail(adjmat_fgesCPDAG_X, adjmat_trueCPDAG))
  shdmat <- matrix(0,nsim, 9)
  # bestk <- matrix(0,nsim, 3)
  
  for(tt in 1:nsim){
    # X_p <- Permute_X(X, seed = tt)
    # saveRDS(X_p, file = paste0(format(today_date, "%Y-%m-%d"), "-vers-", tt, "-permutation",".rds"))
    # XX <- X_p$Xperm
    XX <- X
    rho.est <- rep(1,p)
    # get lambda path ---------------------------------------------------------
    lambda.max <- rep(0, p)
    for(i in 1:(p-1)) lambda.max[i+1] <- norm(2*t(XX[,1:i])%*%(XX[,i+1]*rho.est[i+1]), type = "i")
    lambda.max <- max(lambda.max)
    lambda.path <- lseq(lambda.max/50, lambda.max, lambda.len)
    write.table(lambda.path, file = "lambda_path.txt")
    # some initializaiton -----------------------------------------------------
    BICscores_main <- BICscores_flipflop <- BICscores_bench <- rep(0, lambda.len)
    minrowcor_main <- minrowcor_bench <- minrowcor_flip <- rep(0,lambda.len)
    s0_seq <- rep(0, lambda.len)
    
    for(k in 1:lambda.len){
      set.seed(1)
      lambda <- lambda.path[k]

      # Main algorithm ----------------------------------------------------------
      result_main <- main_iteration_fast(X = XX, vers = k*2, lambda = lambda,
                                         lambda2 = 0.01, max.iter = max.iter,
                                         zeropos = estimands$zeropos)
      s0_seq[k] <- result_main$s0_est
      saveRDS(result_main, file = paste0(format(today_date, "%Y-%m-%d"), "-nsim-", tt,"-vers-",
                                         k, "-lam-", k, "-mainResult",".rds"))
      
      heatmap.2(abs(solve(result_main$theta_est)), dendrogram = "none",  Rowv = F, Colv = F, 
                main= paste0("main_", k), trace = "none")

      # save plots
      png(paste0("NLikelihood", k, "-tt", ".png"))
      plot(result_main$lik_seq, pch = 16, type = "b", cex = 0.7)
      dev.off()
      # # --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
      cor_est <- cor(t(chol(result_main$theta_est)%*%(XX - XX%*%result_main$b_est)))
      minrowcor_main[k] <- sum(abs(cor_est[upper.tri(cor_est)]))
      mle_result <- dag_mle_estimation(X = XX, Bhat = result_main$b_est, Lhat = chol(result_main$theta_est))
      saveRDS(mle_result, file = paste0(format(today_date, "%Y-%m-%d"),  "-nsim-", tt,
                                        "-vers-", k, "-lam-", k,"-mainMLEResult",".rds"))
      BIC_result <- BIC_dag(X = XX, bmle = mle_result$Bmle, omgmle = mle_result$omgmlesq,
                            theta = result_main$theta_est)
      BICscores_main[k] <- BIC_result$BIC
      

      # benchmark ---------------------------------------------------------------
      result_bench <- lassoIdentTheta(X = XX, lambda = lambda, vers = k*2, zeropos = estimands$zeropos,
                                      max.iter = max.iter, lambda2 = 0.01)
      
      saveRDS(result_bench, file = paste0(format(today_date, "%Y-%m-%d"), "-nsim-", tt,"-vers-",
                                         k, "-lam-", k, "-benchResult",".rds"))
      
      # heatmap.2(result_bench$b_est, dendrogram = "none",  Rowv = F, Colv = F, 
      #           main= paste0("bench_", k), trace = "none")
      # 
      # save plots
      png(paste0("NLikelihood_bench", k, "-tt", ".png"))
      plot(result_bench$likeli_seq, pch = 16, type = "b", cex = 0.7)
      dev.off()
      # # --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
      cor_est <- cor(XX - XX%*%result_bench$b_est)
      minrowcor_bench[k] <- sum(abs(cor_est[upper.tri(cor_est)]))
      mle_result_bench <- dag_mle_estimation(X = XX, Bhat = result_bench$b_est, Lhat = diag(n))
      saveRDS(mle_result_bench, file = paste0(format(today_date, "%Y-%m-%d"),  "-nsim-", tt,
                                        "-vers-", k, "-lam-", k,"-benchMLEResult",".rds"))
      BIC_result <- BIC_dag(X = XX, bmle = mle_result_bench$Bmle, 
                            omgmle = mle_result_bench$omgmlesq, theta = diag(n))
      BICscores_bench[k] <- BIC_result$BIC
      #flipflop -----------------------
      # set.seed(10)
      # n <- args$n
      # p <- estimands$realp
      # resultflip <- flipflop(X = XX, n = args$n, p = estimands$realp,
      #                        zeropos_theta = estimands$zeropos,
      #                        max.iter = 20,
      #                        lambda = lambda.path[k],
      #                        lambda2 = 0.01)
      # 
      # psi_hat <- resultflip$psi_est
      # theta_hat <- resultflip$theta_est
      # psi_zeros <- which(abs(psi_hat) < 1e-4, arr.ind = T)
      # theta_zeros <- which(abs(theta_hat) < 1e-4, arr.ind = T)
      # if(any(diff(resultflip$likeli_seq[resultflip$likeli_seq!=0]) > 0))
      #   warning("likelihood not descreasing ! \n")
      # 
      # png(paste0("B_est_flipflop_lik_", k ,".png"))
      # plot(resultflip$likeli_seq[resultflip$likeli_seq!=0], pch = 16, type = "b")
      # lines(resultflip$likeli_seq[resultflip$likeli_seq!=0], type = "o")
      # dev.off()
      # 
      # cat("[INFO] fliplop MLE calculation. \n")
      # resultMLE <- flipflop(X=XX, n = n, p = p,
      #                       zeropos_theta = theta_zeros,
      #                       zeropos_psi = psi_zeros,
      #                       max.iter = 15,
      #                       lambda = 0.0001,
      #                       lambda2 = 0.0001)
      # psi_mle <- resultMLE$psi_est
      # b_est <- resultMLE$b_est
      # theta_mle <- resultMLE$theta_est
      # 
      # heatmap.2(abs(solve(theta_mle)), 
      #           dendrogram = "none", 
      #           Rowv = F, Colv = F, 
      #           main= paste0("flipflop_", k),
      #           trace = "none")
      # 
      # saveRDS(resultMLE, file = paste0(format(today_date, "%Y-%m-%d"),"-lam-",
      #                                  k,"-flipflopMLE",".rds"))
      # BICscores_flipflop[k] <- -n*log(det(psi_mle)) - p*log(det(theta_mle)) +
      #   sum(diag(theta_mle%*%XX%*%psi_mle%*%t(XX))) + log(max(n,p))*sum(abs(psi_mle) > 1e-4)
      # cor_est <- cor(t(chol(resultMLE$theta_est)%*%(XX - XX%*%b_est)))
      # # sd_fliplfop[k] <- sd(cor_est[upper.tri(cor_est)])
      # minrowcor_flip[k] <- sum(abs(cor_est[upper.tri(cor_est)]))
      
      cat("[INFO] lambda ", k, "is finished.\n")
    }
    
    
    png(paste0("BICscores_main", tt ,".png"))
    plot(BICscores_main, pch = 16)
    dev.off()
    png(paste0("minrowcor_main", tt ,".png"))
    plot(minrowcor_main, pch = 16)
    dev.off()
    png(paste0("BICscores_bench", tt ,".png"))
    plot(BICscores_bench, pch = 16)
    dev.off()
    png(paste0("minrowcor_bench", tt ,".png"))
    plot(minrowcor_bench, pch = 16)
    dev.off()
    png(paste0("BICscores_flipflop", tt ,".png"))
    plot(BICscores_flipflop, pch = 16)
    dev.off()
    png(paste0("minrowcor_flip", tt ,".png"))
    plot(minrowcor_flip, pch = 16)
    dev.off()
    # 
    bestk_bic_main <- which.min(BICscores_main)
    bestk_cor_main <- which.min(minrowcor_main)
    bestk_bic_bench <- which.min(BICscores_bench)
    bestk_cor_bench <- which.min(minrowcor_bench)
    bestk_bic_flop <- which.min(BICscores_flipflop)
    bestk_cor_flop <- which.min(minrowcor_flip)
    
    allShdS <- get_shd_for_several_methods(kmainbic = bestk_bic_main,
                                           kbenchbic = bestk_bic_bench,
                                           kmaincor = bestk_cor_main,
                                           kbenchcor = bestk_cor_bench,
                                           today_date = today_date,
                                           tt = tt,
                                           adjmat_trueCPDAG = adjmat_trueCPDAG,
                                           thresh = 0.2)
    
    
    # shdmat[tt,3] <- bestk_bic_main
    # shdmat[tt,4] <- bestk_cor_main
    # shdmat[tt,7] <- bestk_bic_flop
    # shdmat[tt,8] <- bestk_cor_flop
    # cat("[INFO] Best k (bic) is ", bestk_bic_main, "\n")
    # 1) my estimate of B
   
    #2)
    bestresult_main_bic <- readRDS(file = paste0(format(today_date, "%Y-%m-%d"),
                                                 "-nsim-", tt,
                                                 "-vers-", bestk_bic_main, "-lam-", bestk_bic_main,
                                                 "-mainResult",".rds"))
    XXtheta <- chol(bestresult_main_bic$theta_est) %*% XX
    dimnames(XXtheta) <- dimnames(XX)
    fgs_theta <- fges(df = XXtheta, penaltydiscount = 2.0, maxDegree = -1,
                      faithfulnessAssumed = TRUE, verbose = F)
    adjmat_fgesCPDAG_theta <- get_adjmat_from_fges(fgs_theta$edges, p = p,
                                                   varnames = dimnames(X)[[2]])
    shdXLBIC <- compute_SHD_detail(adjmat_fgesCPDAG_theta, adjmat_trueCPDAG)
    shdmat[tt,1] <- shdXLBIC$myshd
    
    

    # # 2)
    bestresult_main_cor <- readRDS(file = paste0(format(today_date, "%Y-%m-%d"),
                                                 "-nsim-", tt,
                                                 "-vers-", bestk_cor_main, "-lam-", bestk_cor_main,
                                                 "-mainResult",".rds"))
    XXtheta <- chol(bestresult_main_cor$theta_est) %*% XX
    dimnames(XXtheta) <- dimnames(XX)
    fgs_theta <- fges(df = XXtheta, penaltydiscount = 2.0, maxDegree = -1,
                      faithfulnessAssumed = TRUE, verbose = F)
    adjmat_fgesCPDAG_theta <- get_adjmat_from_fges(fgs_theta$edges, p = p,
                                                   varnames = dimnames(X)[[2]])
    shdXXLcor <- compute_SHD_detail(adjmat_fgesCPDAG_theta, adjmat_trueCPDAG)
    shdmat[tt,2] <- shdXXLcor$myshd
    
    
    # # flipflop ----------------------------------------------------------------
    #  1)
    bestresult_flop_bic <- readRDS(file = paste0(format(today_date, "%Y-%m-%d"),"-lam-",
                                            bestk_bic_flop, "-flipflopMLE",".rds"))
    
    XXthetaflip_bic <- chol(bestresult_flop_bic$theta_est) %*% XX
    dimnames(XXthetaflip_bic) <- dimnames(XX)
    fgs_theta_flip_bic <- fges(df = XXthetaflip_bic, penaltydiscount = 2.0, maxDegree = -1,
                      faithfulnessAssumed = TRUE, verbose = F)
    adjmat_fgesCPDAG_theta_bic <- get_adjmat_from_fges(fgs_theta_flip_bic$edges, p = p,
                                                   varnames = dimnames(X)[[2]])
    shdXX_flip_bic <- compute_SHD_detail(adjmat_fgesCPDAG_theta_bic, adjmat_trueCPDAG)  
    saveRDS(shdXX_flip_bic, "shdXX_flip_bic.rds")
    shdmat[tt,5] <- shdXX_flip_bic$myshd
    # 2)
    bestresult_cor_flop <- readRDS(file = paste0(format(today_date, "%Y-%m-%d"),"-lam-",
                                                 bestk_cor_flop, "-flipflopMLE",".rds"))
    XXthetaflip_cor <- chol(bestresult_cor_flop$theta_est) %*% XX
    dimnames(XXthetaflip_cor) <- dimnames(XX)
    fgs_theta_flip_cor <- fges(df = XXthetaflip_cor, penaltydiscount = 2.0, maxDegree = -1,
                           faithfulnessAssumed = TRUE, verbose = F)
    adjmat_fgesCPDAG_theta_cor <- get_adjmat_from_fges(fgs_theta_flip_cor$edges, p = p,
                                                   varnames = dimnames(X)[[2]])
    shdXX_flip_cor <- compute_SHD_detail(adjmat_fgesCPDAG_theta_cor, adjmat_trueCPDAG)  
    saveRDS(shdXX_flip_cor,"shdXX_flip_cor.rds")
    shdmat[tt,6] <- shdXX_flip_cor$myshd
    cat("[INFO] This is the ", tt, "th simulation!\n")
  }
  return(list(shdmat = shdmat,
              shdX = shdX))
}




# createComparisonplot ----------------------------------------------------
shdComparePlot <- function(allshd){
  shddata <- data.frame(error_type = factor(rep(names(allshd$shdGES),length(allshd)), levels = c("FN", "wrong_dir", 
                                                       "FP", "ud", "du", "myshd"), ordered = T), 
                        count = c(allshd$shdGES,
                                  allshd$shdXmain,
                                  # allShdS$shdXmainCor,
                                  allshd$shdXbenchCor,
                                  allshd$shdXbench,
                                  allshd$shdXmainCor
                                  # unlist(shdXbench)
                                  # unlist(shdXX_flip_bic),
                                  # unlist(shdXX_flip_cor)
                                  # unlist(shdXbenchCor)
                                  # unlist(shdXX_flip_bic)
                        ),
                        method = rep(c("GES",
                                       "mainbic",
                                       # "GESLcor",
                                       # "GESLbic",
                                       # "GESLFLflipbic",
                                       # "GESLflipcor"),
                                       # "maincor",
                                       "benchcor",
                                       "benchbic",
                                     "benchcor"),
                                     # "GESLflipBic"), 
                                     each = 6))
  saveRDS(shddata, "shddata.rds")
}

shddata <- data.frame(error_type = factor(rep(names(shdX),4), 
                                          levels = c("FN", "wrong_dir", 
                                                     "FP", "ud", "du", 
                                                     "myshd"), 
                                          ordered = T), 
                      count = c(allshdave$shdXave,
                                allshdave$shdXbench,
                                allshdave$shdLbic,
                                allshdave$shdXmain
                                ),
                      method = rep(c("GES",
                                     "bench",
                                     "GESL",
                                     "main"),
                                   each = 6))
saveRDS(shddata, "shddata.rds")

png(filename = "SHD1.png")
ggplot(shddata, aes(x=error_type, y = count, fill = method)) +
  geom_bar(stat = 'identity', position = "dodge") + 
  geom_text(aes(label = count), position = position_dodge(width = 1), vjust=-1)
dev.off()
  
png("EstTheta.png")
heatmap.2(abs(solve(bestresult_main_bic$theta_est)), 
          dendrogram = "none", 
          Rowv = F, Colv = F, 
          trace = "none")
dev.off()
