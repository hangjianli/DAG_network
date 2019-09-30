flipflop_onpath <- function(X, estimands, lambda.path){
  # Run all models once ---------------------------------------------------
  today_date <- Sys.Date()
  n <- dim(X)[1]
  p <- dim(X)[2]
  # BICscores_flip <- rep(0, length(lambda.path))
  # calculate lambda1 path from X----------------------------------------------------------
  for(i in 1:length(lambda.path)){
    resultflip <- flipflop(X = X, n = n, p = p, 
                           zeropos = estimands$zeropos, 
                           max.iter = 50, 
                           lambda = lambda.path[i],
                           lambda2 = 0.2)  
    saveRDS(resultflip, file = paste0(format(today_date, "%Y-%m-%d"),"-lam-", i,"-flipflopResult",".rds"))
    png(paste0("B_est_flipflop_lik_", i ,".png"))
    plot(resultflip$likeli_seq[resultflip$likeli_seq!=0], pch = 16, cex = 0.5)
    dev.off()
    # mle_result <- dag_mle_estimation(X = BIC$X, Bhat = resultflip$b_est, Lhat = chol(resultflip$theta_est))
    # BIC_result <- BIC_dag(X = X, bmle = mle_result$Bmle, omgmle = mle_result$omgmlesq, theta = resultflip$theta_est)
    # BICscores_flip[i] <- BIC_result$BIC
  }
  
  return(list(X=X))
}