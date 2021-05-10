setwd("~/Documents/research/dag_network")
args <- readRDS("~/Documents/research/dag_network/output/121006/args.rds")
estimands <- readRDS("~/Documents/research/dag_network/output/121006/estimands.rds")

for(sim in 1:10){
  # dir.create(path = paste0("output/",args$setting, "/", args$setting, "--", sim))
  setwd(paste0("output/",args$setting, "/", args$setting, "--", sim))
  con <- file("test2.log")
  sink(con, append=TRUE)
  sink(con, append=TRUE, type="message")
  cat("===============\n")
  cat(paste0('starting simulation ', sim, '\n'))
  X_ <- readRDS(paste0("~/Documents/research/dag_network/output/",args$setting, "/", 
                       args$setting, "--", sim, '/X.rds'))
  X <- X_$X
  networkDAG_sol_path2(
    X = X,
    block_size=args$block_size,
    zeropos_list = estimands$zeropos_list,
    lambda_len = 10,
    lambda1_max_div = 100,
    lambdaff_max = 100,
    lambdaff_min = 1,
    maxIter = 30,
    lambda2 = 100
  )
  sink() 
  sink(type="message")
  setwd("~/Documents/research/dag_network")
}
simID = args$setting
process_output_ordered(simID = simID, estimands = estimands, thr = 0.1)
get_all_shd_ordered(simID = simID, estimands, args$num_sim)
get_average_shd_ordered(simID = simID, nsim = as.numeric(args$num_sim))

networkDAG_sol_path2 <- function(
  X,
  block_size,
  zeropos_list,
  block_idx=NULL,
  lambda_len=10,
  lambda2=100,
  lambda1_max_div = 10,
  lambdaff_max = 200,
  lambdaff_min = 1,
  maxIter=100
){
  n <- dim(X)[1]
  p <- dim(X)[2]
  
  lambdaff.path <- lseq(lambdaff_min, lambdaff_max, lambda_len)
  
  
  saveRDS(lambdaff.path, file = "lambdaff.path.rds")
  BICscores_ff <- minrowcor_ff <- rep(0, length(lambdaff.path))
  for(k in 1:length(lambdaff.path)){
    set.seed(1)
    cat(paste0("===============================================================","\n"))
    cat(paste0('[INFO] Lambda k = ', k, "\n"))
    
    resff <- flipflop(
      X = X,
      block_size = block_size,
      zeropos_list = zeropos_list,
      block_idx = block_idx,
      lambda1 = lambdaff.path[k],
      lambda2 = 100, 
      maxIter = maxIter
    )
    
    if(any(apply(resff$bhat, 2, function(x) sum(abs(x) > 1e-4)) >= n)){
      cat("[INFO] lambdaff is too small, algorithm stops at k = ", k, "\n")
    }
    
    saveRDS((resff), file = paste0("ff_lam_", k, ".rds"))
    # compute min correlation
    cor_est_ff <- cor(t(chol(resff$thetahat)%*%(X - X%*%resff$bhat)))
    minrowcor_ff[k] <- sum(abs(cor_est_ff[upper.tri(cor_est_ff)]))
    #compute MLE
    mle_result_ff <- dag_mle_estimation(
      X = X,
      Bhat = resff$bhat,
      Lhat = chol(resff$thetahat)
    )
    saveRDS(mle_result_ff, file = paste0("ffMLE_lam_", k, ".rds"))
    # compute BIC using MLE estimates
    BIC_ff_result <- BIC_dag(
      X = X,
      block_idx = block_idx,
      bmle = mle_result_ff$Bmle,
      omgmle = mle_result_ff$omgmlesq,
      theta = resff$thetahat
    )
    BICscores_ff[k] <- BIC_ff_result$BIC
    saveRDS(BIC_ff_result,  paste0('BIC_ff_result_', k, '.rds'))
  }
  saveRDS(BICscores_ff, "BICscores_ff.rds")
  saveRDS(minrowcor_ff, "minrowcor_ff.rds")
}
