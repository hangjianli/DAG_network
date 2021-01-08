 run_bcd <- function(
  X, 
  block_size,
  zeropos_list,
  block_idx = NULL,
  lambda1=1, 
  lambda2=1,
  tol=1e-7,
  maxIter=100
){
  n <- dim(X)[1]
  p <- dim(X)[2]
  thetahat <- diag(n)
  num_blocks <- length(zeropos_list)
  X_iid <- X[seq(1, n, by = block_size),]
  X_iid <- X
  omega2_hat_iid <- estimate_omega_square(X_iid)
  old_hbeta <- matrix(0, p, p)
  old_htheta <- matrix(0, n, n)
  diff_beta <- 1
  diff_theta <- 1
  iter <- 0
  loss_his <- vector(length = maxIter)
  
  correct_order = NULL
  if (!is.null(block_idx)){
    correct_order <- match(rownames(X), names(unlist(block_idx)))
    # paste0('[INFO] Iter: ', iter, "\n"))
    # cat(rownames(X), names(unlist(block_idx)) )  
  }
  
  while((diff_beta > tol || diff_theta > tol) && iter < maxIter){
    bhat <- estimate_b(
      n = n, p = p, X = X,
      theta_hat = thetahat,
      lambda= rep(lambda1, p)
    )
    S <- get_sample_cov(X, sqrt(omega2_hat_iid), bhat)
    thetahat <- estimate_theta(
      S = S, p = p, 
      lambda2 = lambda2,
      block_size = block_size,
      zeropos_list = zeropos_list,
      block_idx = block_idx,
      correct_order=correct_order
    )
    # err_beta <- norm(estimands$b - bhat, type = '2')^2 / p^2 
    # err_theta <- norm(estimands$theta - thetahat, type = '2')^2 / n^2
    curloss <- eval_loss_func(X, sqrt(omega2_hat_iid), bhat, thetahat, S, lambda1, lambda2)
    cat(paste0('[INFO] Iter: ', iter, "\n"))
    cat(paste0('[INFO] Loss: ', round(curloss, 7), "\n"))
    # cat(paste0('[INFO] err_beta: ', round(err_beta,7), "\n"))
    # cat(paste0('[INFO] err_theta: ', round(err_theta,7), "\n"))
    cat(paste0('[INFO] diff_beta: ', round(diff_beta, 7), "\n"))
    cat(paste0('[INFO] diff_theta: ', round(diff_theta, 7), "\n"))
    cat(paste0("---------------------------------------------","\n"))
    diff_beta <- norm(bhat - old_hbeta, type = 'f')^2 / p^2
    diff_theta <- norm(thetahat - old_htheta, type = 'f')^2 / n^2
    old_hbeta <- bhat
    old_htheta <- thetahat
    iter <- iter + 1
    loss_his[iter] <- curloss
    if(iter == 1){
      cat("[INFO] Saving estimates after one iteration. \n")
      bhat_1iter <- bhat
      thetahat_1iter <- thetahat
    }
    Sys.sleep(0.01)
    flush.console()
  }
  # cat(paste0('[INFO] diff_beta: ', round(diff_beta, 7), "\n"))
  # cat(paste0('[INFO] diff_theta: ', round(diff_theta, 7), "\n"))
  dimnames(bhat) <- list(dimnames(X)[[2]], dimnames(X)[[2]])
  return(list(
    bhat=bhat, 
    thetahat=thetahat,
    bhat_1iter=bhat_1iter,
    thetahat_1iter=thetahat_1iter,
    losses=loss_his[loss_his!=0]
  ))
 }

get_sample_cov <- function(
  X,
  h_omega, 
  h_beta
){
  n <- dim(X)[1]
  p <- dim(X)[2]
  rho_root <- 1 / h_omega
  test.temp <- X%*%(diag(rho_root) - sweep(h_beta, MARGIN = 2, FUN = '*', STATS = rho_root))
  test.res <- apply(test.temp,2,tcrossprod)
  S <- matrix(rowSums(test.res), n, n)
  S.scale <- S/p
  return(S.scale)
}

estimate_b <- function(
  n, p, X, 
  theta_hat=diag(n),
  lambda
){
  choles <- chol(theta_hat) 
  est_B <- matrix(0, p, p) 
  XX <- cbind(0, X)
  B_temp <- foreach(
    i=2:p,
    .combine = "cbind",
    .packages = "glmnet"
    ) %dopar% {
    lars.temp <- glmnet(
      x = as.matrix(choles)%*%XX[,1:i],
      y = as.matrix(choles)%*%XX[,i+1],
      alpha = 1,
      intercept = F,
      lambda = lambda[i]/(2*n),
      standardize = F,
      family = "gaussian",
      thresh = 1e-10
    )
    c(as.numeric(lars.temp$beta)[-1], rep(0, p - i + 1))
    }
  est_B[,2:p] <- B_temp
  return(est_B)
}


estimate_theta <- function(
  S,
  p,
  lambda2, # this should be chosen by grid search
  block_size,
  block_idx=NULL,
  zeropos_list=NULL,
  correct_order=NULL,
  seed=1
){
  n <- dim(S)[1]
  set.seed(seed)
  if(!is.null(zeropos_list)){
    num_blocks = length(zeropos_list)  
  }else{
    num_blocks = length(block_idx)  
  }
  sig.blocks <- vector(mode = "list", length = num_blocks)
  if(!is.null(block_idx)){
    for(i in 1:num_blocks){
      zeros = zeropos_list[[i]]
      if(is.null(zeros) || dim(zeros)[1] == 0)
        zeros = NULL
      cat('[INFO] i ', i)
      temp_sig <- glasso(s = S[block_idx[[i]], block_idx[[i]]],
                         thr = 1.0e-7,
                         rho = lambda2/p,
                         zero = zeros,
                         penalize.diagonal = T)
      sig.blocks[[i]] <- cov2cor(temp_sig$w)
    }
  }else{
    for(i in 1:num_blocks){
      zeros = zeropos_list[[i]]
      if(dim(zeros)[1] == 0)
        zeros = NULL
      temp_sig <- glasso(s = S[(1+block_size*(i-1)) : min(block_size*i,n), 
                               (1+block_size*(i-1)) : min(block_size*i,n)],
                         thr = 1.0e-7,
                         rho = lambda2/p,
                         zero = zeros,
                         penalize.diagonal = T)
      sig.blocks[[i]] <- cov2cor(temp_sig$w)
      # sig.blocks[[i]] <- temp_sig$w
    }    
  }
  sig.est <- as.matrix(do.call(bdiag, sig.blocks))
  if(!is.null(correct_order)){
    sig.est <- sig.est[correct_order, correct_order]    
  }
  theta_est <- solve(sig.est)
  return(theta_est)
}



estimate_omega_square <- function(X){
  p <- dim(X)[2]
  n <- dim(X)[1]
  res = numeric(length = p)
  res[1] = 1
  for(i in 2:p){
    x = as.matrix(cbind(0,X[,1:(i-1)]))
    y = as.matrix(X[,i])
    lambda <- sqrt(log(p)/n)
    model = glmnet(x = x,
                 y = y,
                 alpha = 1,
                 lambda = lambda,
                 intercept = F,
                 standardize = F,
                 family = "gaussian",
                 thresh = 1e-10) 
    res[i] <- norm(y - x%*%model$beta, '2')^2/(2*n) + lambda*norm(model$beta, '1')
  }
  return(res)
}


eval_loss_func <- function(X, homega, hbeta, htheta, S, lam1, lam2){
  n <- dim(X)[1]
  p <- dim(X)[2]
  hphi <- sweep(hbeta, 2, homega^2, '/')
  result <- (
    -p*log(det(htheta)) + 
    p*sum(diag(S%*%htheta)) +
    lam1*sum(abs(hphi)) + 
    lam2*sum(abs(htheta))
  )
  return(result)  
}


networkDAG_sol_path <- function(
  X,
  block_size,
  zeropos_list,
  block_idx=NULL,
  lambda_len=10,
  maxIter=100
){
  n <- dim(X)[1]
  p <- dim(X)[2]
  lambda.path <- rev(get_lam_path(p, X, rho.est = rep(1,p), lambda_len, 100))
  saveRDS(lambda.path, file = "lambda_path.rds")
  BICscores_main <- minrowcor_main <- rep(0, length(lambda.path))
  BICscores_1iter_main <- minrowcor_1iter_main <- rep(0, length(lambda.path))
  for(k in 1:length(lambda.path)){
    set.seed(1)
    cat(paste0("===============================================================","\n"))
    cat(paste0('[INFO] Lambda k = ', k, "\n"))
    res <- run_bcd(
      X = X, 
      block_size = block_size,
      zeropos_list = zeropos_list,
      block_idx = block_idx,
      lambda1 = lambda.path[k],
      lambda2 = .01,
      maxIter = maxIter,
      tol = 1e-7)
    # check if max degree s < n
    if(any(apply(res$bhat, 2, function(x) sum(abs(x) > 1e-4)) >= n)){
      cat("[INFO] lambda is too small, algorithm stops at k = ", k, "\n")
      break
    }
    saveRDS(res, file = paste0("main_lam_", k, ".rds"))
    cor_est <- cor(t(chol(res$thetahat)%*%(X - X%*%res$bhat)))
    minrowcor_main[k] <- sum(abs(cor_est[upper.tri(cor_est)]))
    cor_est_1iter <- cor(t(chol(res$thetahat_1iter)%*%(X - X%*%res$bhat_1iter)))
    minrowcor_1iter_main[k] <- sum(abs(cor_est_1iter[upper.tri(cor_est_1iter)]))
    mle_result <- dag_mle_estimation(
      X = X, 
      Bhat = res$bhat, 
      Lhat = chol(res$thetahat)
    )
    mle_result_1iter <- dag_mle_estimation(
      X = X, 
      Bhat = res$bhat_1iter,
      Lhat = chol(res$thetahat_1iter)
    )
    saveRDS(mle_result, file = paste0("mainMLE_lam_", k, ".rds"))
    saveRDS(mle_result_1iter, file = paste0("mainMLE_1iter_lam_", k, ".rds"))
    BIC_result <- BIC_dag(
      X = X,
      bmle = mle_result$Bmle, 
      omgmle = mle_result$omgmlesq,
      theta = res$thetahat
    )
    BICscores_main[k] <- BIC_result$BIC
    BIC_1iter_result <- BIC_dag(
      X = X,
      bmle = mle_result_1iter$Bmle, 
      omgmle = mle_result_1iter$omgmlesq,
      theta = res$thetahat_1iter
    )
    BICscores_1iter_main[k] <- BIC_1iter_result$BIC
  }
  saveRDS(BICscores_main, "BICscores_main.rds")
  saveRDS(minrowcor_main, "minrowcor_main.rds")
  saveRDS(BICscores_1iter_main, "BICscores_1iter_main.rds")
  saveRDS(minrowcor_1iter_main, "minrowcor_1iter_main.rds")
}

bench_sol_path <- function(
  X,
  block_size,
  zeropos_list=zeropos_list,
  lambda_len=10,
  maxIter=100
){
  n <- dim(X)[1]
  p <- dim(X)[2]
  lambda.path <- rev(get_lam_path(p, X, rho.est = rep(1,p), lambda_len, 100))
  saveRDS(lambda.path, file = "lambda_path.rds")
  BICscores_bench <- minrowcor_bench <- rep(0, length(lambda.path))
  for(k in 1:length(lambda.path)){
    set.seed(1)
    cat(paste0("=======================================================","\n"))
    cat(paste0('[INFO] Lambda k: ', k, "\n"))
    res <- run_bcd(
      X = X, 
      baseline_flag = T,
      block_size = block_size,
      zeropos_list = zeropos_list,
      lambda1 = lambda.path[k],
      lambda2 = .001,
      maxIter = maxIter,
      tol = 1e-7)
    
    # check if max degree s < n
    if(any(apply(res$bhat, 2, function(x) sum(abs(x) > 1e-4)) >= n)){
      cat("[INFO] lambda is too small, algorithm stops at k = ", k, "\n")
      break
    }
    saveRDS(res, file = paste0("baseline_lam_", k, ".rds"))
    cor_est <- cor(t(chol(res$thetahat)%*%(X - X%*%res$bhat)))
    minrowcor_bench[k] <- sum(abs(cor_est[upper.tri(cor_est)]))
    mle_result <- dag_mle_estimation(
      X = X, 
      Bhat = res$bhat, 
      Lhat = chol(res$thetahat)
    )
    saveRDS(mle_result, file = paste0("baselineMLE_lam_", k, ".rds"))
    BIC_result <- BIC_dag(
      X = X,
      bmle = mle_result$Bmle, 
      omgmle = mle_result$omgmlesq,
      theta = res$thetahat
    )
    BICscores_bench[k] <- BIC_result$BIC
    saveRDS(BICscores_bench, "BICscores_baseline.rds")
    saveRDS(minrowcor_bench, "minrowcor_baseline.rds")
  }
}

sim_newalgo_ordered <- function(
  args, 
  estimands,
  start_sim, 
  end_sim,
  lamLen=10
){
  for(sim in start_sim:end_sim){
    dir.create(path = paste0("output/",args$setting, "/", args$setting, "--", sim))
    setwd(paste0("output/",args$setting, "/", args$setting, "--", sim))
    X_ <- sim_X(vers = sim, 
                n = args$n,
                p = estimands$realp,
                omg.sq = estimands$omg.sq,
                sig = estimands$sig, 
                b = estimands$b)
    saveRDS(X_, file = "X.rds")
    X <- X_$X
    networkDAG_sol_path(
      X = X, 
      block_size=args$block_size, 
      zeropos_list = estimands$zeropos_list,
      lambda_len = lamLen
    )
    setwd("~/Documents/research/dag_network")
  }
}


get_shd_ordered <- function(
  kmainbic = NULL, 
  kmaincor = NULL, 
  kbenchbic = NULL, 
  kbenchcor = NULL,
  sim,
  simID,
  bstar_adj, 
  s0, 
  b0,
  theta0,
  thresh = 0.1
){
  output_ordered <- vector(mode = "list")
  # main --------------------------------------------------------------------
  n <- dim(theta0)[1]
  p <- dim(b0)[1]
  # BIC
  main_best_bic <- readRDS(file = paste0(simID, '--', sim, '/main_lam_', kmainbic, '.rds'))
  bhat_adj <- 1*(abs(main_best_bic$bhat) > thresh)
  shdXmain <- compute_SHD_dag(adj1 = bhat_adj, adj_true = bstar_adj, s0) %>% unlist()
  shdXmain['beta_l2err'] = norm(b0 - main_best_bic$bhat, type = '2')^2 / s0
  shdXmain['theta_l2err'] = norm(theta0 - main_best_bic$thetahat, type = '2')^2 / n^2
  output_ordered$shdXmain=shdXmain
  # minCor
  main_best_cor <- readRDS(file = paste0(simID, '--', sim, '/main_lam_', kmaincor, '.rds'))
  bhat_adj <- 1*(abs(main_best_cor$bhat) > thresh)
  shdXmainCor <- compute_SHD_dag(adj1 = bhat_adj, adj_true = bstar_adj, s0) %>% unlist()
  shdXmainCor['beta_l2err'] = norm(b0 - main_best_cor$bhat, type = '2')^2 / s0
  shdXmainCor['theta_l2err'] = norm(theta0 - main_best_cor$thetahat, type = '2')^2 / n^2
  output_ordered$shdXmainCor = shdXmainCor
  # baseline ----------------------------------------------------------------
  # BIC
  baseline_best_bic <- readRDS(file = paste0(simID, '--', sim, '/main_lam_', kbenchbic, '.rds'))
  bhat_adj <- 1*(abs(baseline_best_bic$bhat_1iter) > thresh)
  shdXbaseline <- compute_SHD_dag(adj1 = bhat_adj, adj_true = bstar_adj, s0) %>% unlist()
  shdXbaseline['beta_l2err'] = norm(b0 - baseline_best_bic$bhat, type = '2')^2 / s0
  shdXbaseline['theta_l2err'] = norm(theta0 - baseline_best_bic$thetahat, type = '2')^2 / n^2
  output_ordered$shdXbaseline=shdXbaseline
  # minCor
  baseline_best_cor <- readRDS(file = paste0(simID, '--', sim, '/main_lam_', kbenchcor, '.rds'))
  bhat_adj <- 1*(abs(baseline_best_cor$bhat_1iter) > thresh)
  shdXbaselineCor <- compute_SHD_dag(adj1 = bhat_adj, adj_true = bstar_adj, s0) %>% unlist()
  shdXbaselineCor['beta_l2err'] = norm(b0 - baseline_best_cor$bhat, type = '2')^2 / s0
  shdXbaselineCor['theta_l2err'] = norm(theta0 - baseline_best_cor$thetahat, type = '2')^2 / n^2
  output_ordered$shdXbaselineCor = shdXbaselineCor
  return(output_ordered)  
}


process_output_ordered <- function(
  simID='001', 
  estimands,
  thr=0.1
){
  setwd(paste0('output/', simID))
  args <- readRDS('args.rds')
  estimands <- readRDS('estimands.rds')
  bstar_adj <- 1*(abs(estimands$b) > 0)
  for(sim in 1:args$num_sim){
    BICscores_main <- readRDS(paste0(simID, '--', sim, '/BICscores_main.rds'))  
    minrowcor_main <- readRDS(paste0(simID, '--', sim, '/minrowcor_main.rds'))  
    BICscores_baseline <- readRDS(paste0(simID, '--', sim, '/BICscores_1iter_main.rds'))  
    minrowcor_baseline <- readRDS(paste0(simID, '--', sim, '/minrowcor_1iter_main.rds'))  
    bestk_bic_main <- which.min(BICscores_main)
    bestk_cor_main <- which.min(minrowcor_main)
    bestk_bic_baseline <- which.min(BICscores_baseline)
    bestk_cor_baseline <- which.min(minrowcor_baseline)
    SHD_stats <- get_shd_ordered(
      kmainbic = bestk_bic_main, 
      kmaincor = bestk_cor_main, 
      kbenchbic = bestk_bic_baseline, 
      kbenchcor = bestk_cor_baseline,
      sim = sim,
      simID = simID,
      bstar_adj = bstar_adj, 
      s0 = estimands$s0, 
      theta0 = estimands$theta,
      b0 = estimands$b,
      thresh = thr
    )
    saveRDS(SHD_stats, file = paste0(simID, '--', sim, "/SHDstats.rds"))
  }
  setwd("~/Documents/research/dag_network")
}


get_all_shd_ordered <- function(
  simID, 
  estimands,
  num_sim
){
  thrs <- seq(0, 0.5,length.out = 10)
  allShdS <- vector(mode = "list", length = length(thrs))
  bstar_adj <- 1*(abs(estimands$b) > 0)
  for(sim in 1:num_sim){
    setwd(paste0("output/",simID))
    allshd <- readRDS(paste0(simID, "--", sim, '/SHDstats.rds'))
    BICscores_bench <- readRDS(paste0(simID, "--", sim, "/BICscores_1iter_main.rds"))
    minrowcor_bench <- readRDS(paste0(simID, "--", sim, "/minrowcor_1iter_main.rds"))
    BICscores_main <- readRDS(paste0(simID, "--", sim, "/BICscores_main.rds"))
    minrowcor_main <- readRDS(paste0(simID, "--", sim, "/minrowcor_main.rds"))
    bestk_bic_main <- which.min(BICscores_main)
    bestk_cor_main <- which.min(minrowcor_main)
    bestk_bic_baseline <- which.min(BICscores_bench)
    bestk_cor_baseline <- which.min(minrowcor_bench)
    for(j in 1:length(thrs)){
      allShdS[[j]] <- get_shd_ordered(
        kmainbic = bestk_bic_main, 
        kmaincor = bestk_cor_main, 
        kbenchbic = bestk_bic_baseline, 
        kbenchcor = bestk_cor_baseline,
        sim = sim,
        simID = simID,
        bstar_adj = bstar_adj, 
        s0 = estimands$s0, 
        theta0 = estimands$theta,
        b0 = estimands$b,
        thresh = thrs[j]
      )
    }
    bestbenchBIC <- which.min(sapply(
      allShdS, 
      function(x) abs(x$shdXbaseline['pnum'] - allshd$shdXmain['pnum'])
    ))
    bestbenchCor <- which.min(sapply(
      allShdS, 
      function(x) abs(x$shdXbaselineCor['pnum'] - allshd$shdXmainCor['pnum'])
    ))
    allshd$shdXbaseline <- allShdS[[bestbenchBIC]]$shdXbaseline
    allshd$shdXbaselineCor <- allShdS[[bestbenchCor]]$shdXbaselineCor
    saveRDS(allshd, paste0(simID, "--", sim, "/SHDclose.rds"))
    cat("[INFO] Sim", sim, "is done. \n")
    setwd("~/Documents/research/dag_network")
  }
}

# get_all_shd_unordered <- function(
#   simID, 
#   estimands,
#   num_sim
# ){
#   return()
# }


get_average_shd_ordered <- function(simID, nsim){
  SHDres <- readRDS(paste0("output/", simID, '/', simID, '--1/', "SHDclose.rds"))
  num_statistic <- length(names(SHDres[[1]]))
  total <- data.frame(
    row.names = names(SHDres[[1]]), 
    shdXmain=rep(0, num_statistic),
    shdXbaseline=rep(0, num_statistic),
    shdXmainCor=rep(0, num_statistic),
    shdXbaselineCor=rep(0, num_statistic)
  )
  for (sim in 1:nsim){
    allshd <- as.data.frame(
      readRDS(paste0("output/", simID,  '/', simID, "--", sim, "/SHDclose.rds"))
    )
    allshd <- allshd[, names(total)]
    total <- total + allshd
  }
  total <- round(total / nsim, 7)
  total$B0 = c(estimands$s0, rep(0,num_statistic-1))
  saveRDS(total, paste0("output/", simID, "/shd_average.rds"))
}

get_average_shd_unordered <- function(simID, nsim){
  SHDres <- readRDS(paste0("output/", simID, '/', simID, '--1/', "SHDstats.rds"))
  num_statistic <- length(names(SHDres[[1]]))
  total <- data.frame(
    row.names = names(SHDres[[1]]), 
    shdXmain = rep(0, num_statistic),
    shd_pc = rep(0, num_statistic),
    shd_pc_decor=rep(0, num_statistic),
    shd_ges=rep(0, num_statistic),
    shd_ges_decor=rep(0, num_statistic),
    shd_sbn=rep(0, num_statistic),
    shd_sbn_decor=rep(0, num_statistic)
  )
  for (sim in 1:nsim){
    allshd <- as.data.frame(
      readRDS(paste0("output/", simID,  '/', simID, "--", sim, "/SHDstats.rds"))
    )
    allshd <- allshd[, names(total)]
    total <- total + allshd
  }
  total <- round(total / nsim, 2)
  total$B0 = c(estimands$s0, rep(0,num_statistic-1))
  saveRDS(total, paste0("output/", simID, "/shd_average.rds"))
}


GES_sol <- function(X, decor=F){
  n <- dim(X)[1]
  p <- dim(X)[2]
  set.seed(100)
  fgsX <- tetradrunner(
    algoId = 'fges',
    df = X,
    scoreId = 'sem-bic',
    dataType = 'continuous',
    faithfulnessAssumed=TRUE,
    maxDegree=-1,
    verbose=TRUE
  )
  adjmat_fgesCPDAG_X <- get_adjmat_from_fges(
    fgsX$edges,
    p = p, 
    varnames = dimnames(X)[[2]]
  )
  if(decor){
    saveRDS(adjmat_fgesCPDAG_X, file = 'adjmat_fges_CPDAG_decor.rds')  
  }else{
    saveRDS(adjmat_fgesCPDAG_X, file = 'adjmat_fges_CPDAG.rds')  
  }
}

pc_sol <- function(X, decor=F){
  n <- dim(X)[1]
  p <- dim(X)[2]
  suffstat <- list(C = cor(X), n = n)
  res_pc <- pcalg::pc(
    suffStat = suffstat, 
    indepTest = gaussCItest, 
    alpha = 0.05,
    m.max = 5,
    labels = dimnames(X)[[2]]
  )
  adjmat_pc_CPDAG <- as(res_pc, "amat")
  if(decor){
    saveRDS(adjmat_pc_CPDAG, file = 'adjmat_pc_CPDAG_decor.rds')
  }
  else{
    saveRDS(adjmat_pc_CPDAG, file = 'adjmat_pc_CPDAG.rds')  
  }
}

get_Xdecor <- function(Xp){
  bic_score <- readRDS('BICscores_main.rds')
  bic_score_1iter <- readRDS('BICscores_1iter_main.rds')
  best_bic <- which.min(bic_score)
  best_bic_1iter <- which.min(bic_score_1iter)
  cat('[INFO] lambda index for the best BIC: ',
      best_bic, "|", paste0(best_bic_1iter, " (1iter)"),"\n")
  best_res <- readRDS(paste0('main_lam_', best_bic, '.rds'))
  best_res_1iter <- readRDS(paste0('main_lam_', best_bic_1iter, '.rds'))
  X_decor <- chol(best_res$thetahat) %*% Xp
  X_decor_1iter <- chol(best_res_1iter$thetahat) %*% Xp
  dimnames(X_decor) <- dimnames(Xp)
  dimnames(X_decor_1iter) <- dimnames(Xp)
  return(list(X_decor=X_decor, X_decor_1iter=X_decor_1iter))
}


sparsebn_sol <- function(X, decor=F){
  sbX <- sparsebnData(X, type = "continuous")
  sb_path <- estimate.dag(sbX)
  sol_idx <- select.parameter(sb_path, sbX)
  sol_base <- select(sb_path, index = sol_idx)
  adjmat_sbbase <- get.adjacency.matrix(sol_base) %>% as.matrix()
  if(decor){
    saveRDS(adjmat_sbbase, file = 'adjmat_sparsebn_CPDAG_decor.rds')
  }
  else{
    saveRDS(adjmat_sbbase, file = 'adjmat_sparsebn_CPDAG.rds')  
  }
}


sim_newalgo_unordered <- function(
  args, 
  estimands, 
  start_sim=1, 
  end_sim=args$num_sim, 
  lamLen=15){
  for(sim in start_sim:end_sim){
    dir.create(path = paste0("output/",args$setting, "/", args$setting, "--", sim))
    setwd(paste0("output/",args$setting, "/", args$setting, "--", sim))
    X_ <- sim_X(vers = sim, 
                n = args$n,
                p = estimands$realp,
                omg.sq = estimands$omg.sq,
                sig = estimands$sig, 
                b = estimands$b)
    X <- X_$X
    XXp <- Permute_X(X, seed = sim)
    Xp <- XXp$Xperm
    saveRDS(X_, file = "X.rds")
    saveRDS(XXp, file = "Xp.rds")
    # Other methods -----------------------------------------------------------
    cat('[INFO] Running GES... \n')
    GES_sol(Xp, decor = F)
    cat('[INFO] Running PC... \n')
    pc_sol(Xp, decor = F)
    cat('[INFO] Running sparsebn... \n')
    sparsebn_sol(Xp, decor = F)
    # NetworkDAG --------------------------------------------------------------
    networkDAG_sol_path(
      X = Xp, 
      block_size=args$block_size, 
      zeropos_list = estimands$zeropos_list,
      lambda_len = lamLen,
      maxIter = 100
    )
    # decorrelation -----------------------------------------------------------
    Xdecor <- get_Xdecor(Xp)
    GES_sol(Xdecor, decor = T)
    pc_sol(Xdecor, decor = T)
    sparsebn_sol(Xdecor, decor = T)
    setwd("~/Documents/research/dag_network")  
  }
}


get_shd_unordered <- function(
  bestk_bic_main, 
  sim,
  simID,
  bstar_adj, 
  s0, 
  thresh = 0.1
){
  output_unordered <- vector(mode = "list")
  main_best_bic <- readRDS(file = paste0(simID, '--', sim, '/main_lam_', bestk_bic_main, '.rds'))
  bhat_adj_cpdag <- bnstruct::dag.to.cpdag(1*(abs(main_best_bic$bhat) > thresh))
  shdXmain <- unlist(compute_SHD_detail(bhat_adj_cpdag, bstar_adj, s0))
  # get PC shd --------------------------------------------------------------
  adjmat_pc_CPDAG <- readRDS(
    file = paste0(simID, '--', sim, '/adjmat_pc_CPDAG.rds'))
  adjmat_fges_CPDAG_decor <- readRDS(
    file = paste0(simID, '--', sim, '/adjmat_pc_CPDAG_decor.rds'))
  shd_pc <- unlist(compute_SHD_detail(adjmat_pc_CPDAG, bstar_adj, s0))
  shd_pc_decor <- unlist(compute_SHD_detail(adjmat_fges_CPDAG_decor, bstar_adj, s0))
  # get GES shd -------------------------------------------------------------
  adjmat_fges_CPDAG <- readRDS(
    file = paste0(simID, '--', sim, '/adjmat_fges_CPDAG.rds'))
  adjmat_fges_CPDAG_decor <- readRDS(
    file = paste0(simID, '--', sim, '/adjmat_fges_CPDAG_decor.rds'))
  shd_ges <- unlist(compute_SHD_detail(adjmat_fges_CPDAG, bstar_adj, s0))
  shd_ges_decor <- unlist(compute_SHD_detail(adjmat_fges_CPDAG_decor, bstar_adj, s0))
  # get sparsebn shd --------------------------------------------------------
  adjmat_sparsebn_CPDAG <- readRDS(
    file = paste0(simID, '--', sim, '/adjmat_sparsebn_CPDAG.rds'))
  adjmat_sparsebn_CPDAG_decor <- readRDS(
    file = paste0(simID, '--', sim, '/adjmat_sparsebn_CPDAG_decor.rds'))
  shd_sbn <- unlist(compute_SHD_detail(adjmat_sparsebn_CPDAG, bstar_adj, s0))
  shd_sbn_decor <- unlist(compute_SHD_detail(adjmat_sparsebn_CPDAG_decor, bstar_adj, s0))
  
  output_unordered$shdXmain=shdXmain
  output_unordered$shd_pc=shd_pc
  output_unordered$shd_pc_decor=shd_pc_decor
  output_unordered$shd_ges=shd_ges
  output_unordered$shd_ges_decor=shd_ges_decor
  output_unordered$shd_sbn=shd_sbn
  output_unordered$shd_sbn_decor=shd_sbn_decor
  return(output_unordered)
}


process_output_unordered <- function(simID = simID, nsim, thr = 0.1){
  setwd(paste0('output/', simID))
  estimands <- readRDS('estimands.rds')
  bstar_adj_cpdag <- bnstruct::dag.to.cpdag(1*(estimands$b != 0)) 
  for(sim in 1:nsim){
    cat(paste0('[INFO] Processing sim ', sim, '\n'))
    BICscores_main <- readRDS(paste0(simID, '--', sim, '/BICscores_main.rds'))  
    bestk_bic_main <- which.min(BICscores_main)
    SHD_stats <- get_shd_unordered(
      bestk_bic_main = bestk_bic_main, 
      sim = sim,
      simID = simID,
      bstar_adj = bstar_adj_cpdag, 
      s0 = estimands$s0, 
      thresh = thr
    )
    saveRDS(SHD_stats, file = paste0(simID, '--', sim, "/SHDstats.rds"))
  }
  setwd("~/Documents/research/dag_network")
}

sample_sc_data <- function(
  full_log_vals,
  full_idx,
  size=rep(20, 7),
  seed=1
){
  set.seed(seed)  
  subsetID <- sort(unlist(mapply(sample, full_idx, size=size)))
  subsetXp <- full_log_vals[unlist(subsetID), ]
  temp <- substr(dimnames(subsetXp)[[1]], start = 0, stop = 3) 
  cellnames <- sub(pattern = "_$",replacement =  "", x = temp)
  return(list(
    subsetXp=subsetXp,
    cellnames=cellnames
  ))
}

get_cell_block_idx <- function(cellnames){
  #' Return a list where each element is a vector of row indices belonging to the same block/cluster
  #' 
  unique_cell <- unique(cellnames)
  cellname_idx <- numeric(length = length(unique_cell))
  for(i in 1:length(unique_cell)){
    cellname_idx[i] <- which(cellnames == unique_cell[i]) %>% head(1)
  }
  cellname_idx <- c(cellname_idx, length(cellnames)+1)
  block_idx <- list()
  for(i in 1:length(unique_cell)){
    block_idx[[i]] <- seq(cellname_idx[i], cellname_idx[i+1]-1)
  }
  return(block_idx)
}

plot_cpdag <- function(cpdag){
  idx <- which(cpdag !=0, arr.ind=T)
  edgelist <- cbind(rownames(cpdag)[idx[,"row"]], colnames(cpdag)[idx[,"col"]])
  g <- graph_from_edgelist(edgelist)
  plot(g, edge.arrow.size=0.5, vertex.size=5, 
       vertex.label.dist=1.5,
       label.cex=.1)
}

