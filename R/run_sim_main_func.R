run_sim_main_func <- function(args, max.iter = 30, seed = 10, lambda2 = 0.2,
                         lambda.len = 10, threshold = seq(0,0.2,length.out = 4)){
  today_date <- Sys.Date()
  estimands <- generate_parameters(args, seed = seed)
  png("B_true.png")
  # 2. Create the plot
  heatmap.2(abs(estimands$b), dendrogram = "none", Rowv = F, Colv = F, trace = "none")
  # 3. Close the file
  dev.off()
  saveRDS(estimands, file = paste0(format(today_date, "%Y-%m-%d"), "-vers-", args$setting, "-estimands",".rds"))
  saveRDS(args, file = paste0(format(today_date, "%Y-%m-%d"), "-vers-", args$setting, "-args",".rds"))
  p = estimands$realp
  summary <- vector("list", length = args$num_sim)
  for(w in 1:args$num_sim){
    sim_index = paste0(today_date, '-', args$setting, "-", w)
    # # If X is given, comment out the section below
    X_ <- sim_X(vers = w, p = p,
                args = args, omg.sq = estimands$omg.sq,
                sig = estimands$sig, b = estimands$b)
    X <- X_$X
    # Run all models once ---------------------------------------------------
    n <- dim(X)[1]
    p <- dim(X)[2]
    # calculate lambda1 path from X----------------------------------------------------------
    rho.est <- rep(1,p)
    lambda.max <- rep(0, p)
    for(i in 1:(p-1)) lambda.max[i+1] <- norm(2*t(X[,1:i])%*%(X[,i+1]*rho.est[i+1]), type = "i")
    lambda.max <- max(lambda.max)
    # The starting point is important. Cannot make it to small.
    lambda.path <- lseq(lambda.max/30, lambda.max/2, lambda.len)
    theta_errors <- matrix(0, ncol = 3, nrow = lambda.len) 
    BICscores_main <- rep(0, lambda.len)
    for(k in 1:lambda.len){
      start <- Sys.time()
      lambda <- lambda.path[k]
      cat("Lambda position:", k, "Running main model...", "\n")
      result2 <- main_iteration_fast(X = X, vers = w*k, lambda = lambda, lambda2 = lambda2,
                                     max.iter = max.iter*2, zeropos = estimands$zeropos)
      png(paste0("B_est_main_lik_", k ,".png"))
      plot(result2$lik_seq, pch = 16, cex = 0.5)
      dev.off()
      mle_result <- dag_mle_estimation(X = X, Bhat = result2$b_est, Lhat = chol(result2$theta_est))
      BIC_result <- BIC_dag(X = X, bmle = mle_result$Bmle, omgmle = mle_result$omgmlesq, theta = result2$theta_est)
      BICscores_main[k] <- BIC_result$BIC
      cat("[INFO] lambda index ", k, " BIC score: ", BICscores_main[k], "\n")
      saveRDS(result2, file = paste0(format(today_date, "%Y-%m-%d"), "-vers-", w, "-lam-", k,"-mainResult",".rds"))
      if(sum(abs(result2$b_est)) > 0){
        png(paste0("B_est_main_", k ,".png"))
        heatmap.2(abs(result2$b_est), dendrogram = "none", Rowv = F, Colv = F, trace = "none")
        dev.off()
      }
      print("Second model done. \n")
    }
    cat("This is the ", w, "th sample!\n")
  }
  return(list(X=X,
              estimands=estimands,
              BICscores = BICscores_main,
              lambdas = lambda.path))
}