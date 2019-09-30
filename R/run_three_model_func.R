run_three_model <- function(X, vers, estimands, 
                            lambda2 = 0.1, 
                            max.iter = 30,
                            lambda.len = 10, 
                            threshold = seq(0,0.2,length.out = 4)){

  n <- dim(X)[1]
  p <- dim(X)[2]
  today_date <- Sys.Date()
  # calculate lambda1 path from X----------------------------------------------------------
  rho.est <- rep(1,p)
  lambda.max <- rep(0, p)
  for(i in 1:(p-1)) lambda.max[i+1] <- norm(2*t(X[,1:i])%*%(X[,i+1]*rho.est[i+1]), type = "i") 
  lambda.max <- max(lambda.max)
  # The starting point is important. Cannot make it to small.
  lambda.path <- lseq(lambda.max/20, lambda.max/2, lambda.len)
  
  theta_errors <- matrix(0, ncol = 3, nrow = lambda.len) 
  colnames(theta_errors) <- c("lassoIdent", "main", "scaledLasso")
  
  bench_TPFP <- main_TPFP <- twoStep_TPFP <- vector("list", lambda.len)
  berror_bench <- berror_main <- berror_twoStep <- matrix(0, nrow = length(threshold), ncol = lambda.len)
  bench.hamming <- main.hamming <- twoStep.hamming <- vector("list", length = lambda.len)
  bench.edgenum <- main.edgenum <- twostep.edgenum <- matrix(0, nrow = length(threshold), ncol = lambda.len)
  BICscores_main <- rep(0, lambda.len)
  
  for(k in 1:lambda.len){
    start <- Sys.time()
    lambda <- lambda.path[k]
    cat("Lambda position: ", k, ". ", "Running lassoIdentTheta model...", "\n")
    result1 <- lassoIdentTheta(X = X, vers = vers*k, lambda = lambda, lambda2 = lambda2,
                               max.iter = max.iter, zeropos = estimands$zeropos)
    if(sum(abs(result1$b_est)) > 0){
      png(paste0("B_est_bench_", k ,".png"))
      heatmap.2(abs(result1$b_est), dendrogram = "none", Rowv = F, Colv = F, trace = "none")
      dev.off()
    }
    print("First model done. \n")
    # benchmark ----------------------------------------------------------------
    theta_errors[k,"lassoIdent"] <- norm((estimands$theta - result1$theta.est), "f")/(n^2-dim(estimands$zeropos)[1])
    num_FP_bench <- sapply(threshold, function(x) sum(abs(result1$b_est[upper.tri(result1$b_est)][abs(estimands$b[upper.tri(estimands$b)]) < 1e-7]) > x))
    num_TP_bench <- sapply(threshold, function(x) sum(abs(result1$b_est[upper.tri(result1$b_est)][abs(estimands$b[upper.tri(estimands$b)]) > 1e-7]) > x))
    bench.hamming[[k]] <- sapply(threshold, function(x) sum(abs(result1$b_est[upper.tri(result1$b_est)][abs(estimands$b[upper.tri(estimands$b)]) > 1e-7]) <= x)) + num_FP_bench # FN + FP
    bench_TPFP[[k]] <- list(FP = num_FP_bench, TP = num_TP_bench)
    bench.edgenum[,k] <- unlist(lapply(lapply(threshold, function(x) ifelse(abs(result1$b_est) > x, sign(result1$b_est)*(abs(result1$b_est) - x), 0)),
                                       function(x) sum(abs(x) > 1e-7)))
    bench_listofB <- lapply(threshold, function(x) ifelse(abs(result1$b_est) > x, sign(result1$b_est)*(abs(result1$b_est) - x), 0))
    berror_bench[,k] <- sapply(bench_listofB, function(x) norm((x - estimands$b), type = "f") / sqrt(estimands$s0))
    # main --------------------------------------------------------------------
    cat("Lambda position:", k, "Running main model...", "\n")
    result2 <- main_iteration_fast(X = X, vers = vers*k, lambda = lambda, lambda2 = lambda2,
                                   max.iter = max.iter*2, zeropos = estimands$zeropos)
    mle_result <- dag_mle_estimation(X = X, Bhat = result2$b_est, Lhat = chol(result2$theta_est))
    BIC_result <- BIC_dag(X = X, bmle = mle_result$Bmle, omgmle = mle_result$omgmlesq, theta = result2$theta_est)
    BICscores_main[k] <- BIC_result$BIC
    cat("[INFO] lambda index ", k, " BIC score: ", BICscores_main[k], "\n")
    saveRDS(result2, file = paste0(format(today_date, "%Y-%m-%d"), "-vers-", vers, "-lam-", k,"-mainResult",".rds"))
    if(sum(abs(result2$b_est)) > 0){
      png(paste0("B_est_main_", k ,".png"))
      heatmap.2(abs(result2$b_est), dendrogram = "none", Rowv = F, Colv = F, trace = "none")
      dev.off()
    }
    print("Second model done. \n")
    num_FP_main <- sapply(threshold, function(x) sum(abs(result2$b_est[upper.tri(result2$b_est)][abs(estimands$b[upper.tri(estimands$b)]) < 1e-7]) > x))
    num_TP_main <- sapply(threshold, function(x) sum(abs(result2$b_est[upper.tri(result2$b_est)][abs(estimands$b[upper.tri(estimands$b)]) > 1e-7]) > x))
    main_TPFP[[k]] <- list(FP=num_FP_main, TP=num_TP_main)
    main_listofB <- lapply(threshold, function(x) ifelse(abs(result2$b_est) > x, sign(result2$b_est)*(abs(result2$b_est) - x), 0))
    berror_main[,k] <- sapply(main_listofB, function(x) norm((x - estimands$b), type = "f") / sqrt(estimands$s0))
    theta_errors[k,"main"] <- norm((estimands$theta - result2$theta_est), type = "f")/(n^2-dim(estimands$zeropos)[1])
    main.edgenum[,k] <- unlist(lapply(lapply(threshold, function(x) ifelse(abs(result2$b_est) > x, 
                                                                           sign(result2$b_est)*(abs(result2$b_est) - x), 0)), function(x) sum(abs(x) > 1e-7)))
    main.hamming[[k]] <- sapply(threshold, function(x) sum(abs(result2$b_est[upper.tri(result2$b_est)][abs(estimands$b[upper.tri(estimands$b)]) > 1e-7]) <= x)) + num_FP_main
    # two.step ----------------------------------------------------------------
    cat("Lambda position: ", k, ". ", "Running twostep model...", "\n")
    result3 <- scale_lasso_init(X = X, vers = vers*k, lambda = lambda, lambda2 = lambda2, zeropos = estimands$zeropos)
    if(sum(abs(result3$b_est)) > 0){
      png(paste0("B_est_twostep_", k ,".png"))  
      heatmap.2(abs(result3$b_est), dendrogram = "none", Rowv = F, Colv = F, trace = "none")
      dev.off()  
    }
    cat("Two-step model done. \n")
    num_TP_two <- sapply(threshold, function(x) sum(abs(result3$b_est[upper.tri(result3$b_est)][abs(estimands$b[upper.tri(estimands$b)]) > 1e-7]) > x))
    num_FP_two <- sapply(threshold, function(x) sum(abs(result3$b_est[upper.tri(result3$b_est)][abs(estimands$b[upper.tri(estimands$b)]) < 1e-7]) > x))
    twoStep_TPFP[[k]] <- list(FP=num_FP_two, TP=num_TP_two)
    twostep_listofB <- lapply(threshold, function(x) ifelse(abs(result3$b_est) > x, sign(result3$b_est)*(abs(result3$b_est) - x), 0))
    berror_twoStep[,k] <- sapply(twostep_listofB, function(x) norm((x - estimands$b), type = "f") / sqrt(estimands$s0))
    twostep.edgenum[, k] <- unlist(lapply(lapply(threshold, function(x) ifelse(abs(result3$b_est) > x, sign(result3$b_est)*(abs(result3$b_est) -x), 0)),
                                          function(x) sum(abs(x) > 1e-7)))
    theta_errors[k, "scaledLasso"] <- norm((estimands$theta - result3$theta.est), "f")/(n^2-dim(estimands$zeropos)[1])
    twoStep.hamming[[k]] <- sapply(threshold, function(x) sum(abs(result3$b_est[upper.tri(result3$b_est)][abs(estimands$b[upper.tri(estimands$b)]) > 1e-7]) <= x)) + num_FP_two
    
    end <- Sys.time()
    cat("Round ", k, " of the three finished. ")
    print(end-start)
    cat("\n")
  }
  
  return(list(X = X,
              lambda.path = lambda.path,
              threshold = threshold,
              lambda2 = lambda2,
              theta_errors = theta_errors,
              bench_TPFP=bench_TPFP, 
              main_TPFP=main_TPFP,
              twoStep_TPFP=twoStep_TPFP,
              berror_bench=berror_bench,
              berror_main=berror_main,
              berror_twoStep=berror_twoStep,
              bench.hamming=bench.hamming,
              main.hamming=main.hamming,
              twoStep.hamming=twoStep.hamming,
              bench.edgenum=bench.edgenum,
              main.edgenum=main.edgenum,
              twostep.edgenum=twostep.edgenum,
              BICscores_main = BICscores_main))
}

