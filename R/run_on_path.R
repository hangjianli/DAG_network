for(k in 1:length(lambda.path)){
  setwd("~/OneDrive/Documents/research/code")
  
  lambda <- rep(lambda.path[k], p)
  # lambda
  source("lassoIdentThetaIter.R") # without considering theta -- ignoring all row-wise correlation
  # source("benchmark_elas.R")
  if(!all(diff(likeli[likeli!=0]) <= 0)) cat("The likelihood for benchmark is not descreasing !!!")
  # check convergence
  plot(likeli, pch = 16, main = "Likelihood value for benchmark iterations.")
  plot(rho_diff_, pch = 16, main = "Difference in \rho for benchmark.")
  plot(phi_diff_, pch = 16)
  cat("Benchmark done. \n")
# original ----------------------------------------------------------------
  theta_errors[k,"lassoIdent"] <- Theta_error_larsvar
  bench_TPFP[[k]] <- list(num_FP_bench, num_TP_bench)
  # theta_error_lasso[k] <- Theta_error_larsvar
  berror_bench[k] <- B_error_bench
  bench.hamming[[k]] <- hamming.dist2
  bench.edgenum[,k] <- unlist(lapply(lapply(threshold,
                                            function(x) ifelse(abs(b_est2) > x, 
                                                               sign(b_est2)*(abs(b_est2) - x), 0)),
                                     function(x) sum(abs(x) > 1e-7)))
  
# elastic net -------------------------------------------------------------
  # bench_TPFP[[k]] <- list(num_FP_bench_elas, num_TP_bench_elas)
  # berror_bench[k] <- B_error_bench_elas
  # bench.hamming[[k]] <- hamming.dist2_elas
  # bench.edgenum[,k] <- unlist(lapply(lapply(threshold,
  #                                           function(x) ifelse(abs(b_est_elas) > x,
  #                                                              sign(b_est_elas)*(abs(b_est_elas) - x), 0)),
  #                                    function(x) sum(abs(x) > 1e-7)))
  
  # setwd("~/OneDrive/Documents/research/code")
  setwd("C:/Users/hangjian/OneDrive/Documents/research/code")
  
  source("main_iteration_fast.R")
  plot(total.likeli[total.likeli!=0], pch = 16, cex = 0.7, main = paste0("Main", lambda[1]))
  plot(rho_difff[1:iter2-1], pch = 16, cex = 0.7, main = paste0("Main", lambda[1]))
  plot(phi_difff[1:iter2-1], pch = 16, cex = 0.7, main = paste0("Main", lambda[1]))
  plot(theta_difff_[1:iter2-1], pch = 16, cex = 0.7, main = paste0("Main", lambda[1]))
  
  # source("elasticnet.R")
  if(!all(diff(total.likeli[total.likeli!=0]) <= 0)) {
    cat("The likelihood for main is not descreasing !!!", "\n")
    cat("k = ", k, "\n")
    plot(total.likeli, pch = 16, cex = 0.6)
  }
  
# original main -----------------------------------------------------------
  main_TPFP[[k]] <- list(num_FP_main, num_TP_main)
  # main_TPFP[[k]] <- list(num_FP_elas, num_TP_elas) 
  berror_main[k] <- B_error[iter2-1]
  theta_errors[k,"main"] <- Theta_error[iter2-1] 
  main.edgenum[,k] <- unlist(lapply(lapply(threshold,
                                           function(x) ifelse(abs(b_est) > x, sign(b_est)*(abs(b_est) - x), 0)), 
                                    function(x) sum(abs(x) > 1e-7)))
  main.hamming[[k]] <- hamming.dist 
  cat("Theta_error Main: ", Theta_error[iter2 - 1], "\n")
  cat("Main done. \n")
# two.step ----------------------------------------------------------------
  # setwd("~/OneDrive/Documents/research/code")
  setwd("C:/Users/hangjian/OneDrive/Documents/research/code")
  
  source("scaled_lasso_init.R")
  twoStep_TPFP[[k]] <- list(num_FP_two, num_TP_two)
  berror_twpStep[k] <- B_error_heuristic
  twostep.edgenum[, k] <- unlist(lapply(lapply(threshold, function(x) ifelse(abs(b_est_heur) > x, sign(b_est_heur)*(abs(b_est_heur) -x), 0)),
                                        function(x) sum(abs(x) > 1e-7)))
  theta_errors[k, "scaledLasso"] <- Theta_error_heur
  twoStep.hamming[[k]] <- hamming.dist3
  
  cat("Two-step done. \n")
  cat("Lambda: ", k, "\n")
}
   