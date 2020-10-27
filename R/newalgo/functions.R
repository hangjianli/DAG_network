run_bcd <- function(
  #' 
  #' @block_size: args$block_size
  #' @zeropos_list: estimands$zeropos_list
  #' 
  X, 
  block_size,
  estimands,
  zeropos_list,
  baseline_flag=FALSE,
  lambda1=1, 
  lambda2=1,
  tol=1e-7,
  maxIter=100
){
  n <- dim(X)[1]
  p <- dim(X)[2]
  thetahat <- diag(n)
  num_blocks <- ceiling(n/block_size)
  X_iid <- X[seq(1, n, by = block_size),]
  X_iid <- X
  omega2_hat_iid <- estimate_omega_square(X_iid)
  
  old_hbeta <- matrix(0, p, p)
  old_htheta <- matrix(0, n, n)
  diff_beta <- 1
  diff_theta <- 1
  iter <- 0
  while((diff_beta > tol || diff_theta > tol) && iter < maxIter){
    bhat <- estimate_b(
      n = n, p = p, X = X,
      theta_hat = thetahat,
      lambda= rep(lambda1, p)
    )
    S <- get_sample_cov(X, sqrt(omega2_hat_iid), bhat)
    if(baseline_flag){
      cat(paste0("[INFO] Assume no row-wise correlations. \n"))
      return(list(bhat=bhat, thetahat=thetahat))
    }
    thetahat <- estimate_theta(
      S = S, p = p, lambda2 = lambda2,
      num_blocks = num_blocks,
      block_size = block_size,
      zeropos_list = zeropos_list
    )
    err_beta <- norm(estimands$b - bhat, type = '2') 
    err_theta <- norm(estimands$theta - thetahat, type = '2')
    curloss <- eval_loss_func(X, sqrt(omega2_hat_iid), bhat, thetahat, S, lambda1, lambda2)
    cat(paste0('[INFO] Iter: ', iter, "\n"))
    cat(paste0('[INFO] Loss: ', round(curloss, 7), "\n"))
    cat(paste0('[INFO] err_beta: ', round(err_beta,7), "\n"))
    cat(paste0('[INFO] err_theta: ', round(err_theta,7), "\n"))
    cat(paste0('[INFO] diff_beta: ', round(diff_beta, 7), "\n"))
    cat(paste0('[INFO] diff_theta: ', round(diff_theta, 7), "\n"))
    cat(paste0("---------------------------------------------","\n"))
    diff_beta <- norm(bhat - old_hbeta, type = 'f') / p^2
    diff_theta <- norm(thetahat - old_htheta, type = 'f') / n^2
    old_hbeta <- bhat
    old_htheta <- thetahat
    iter <- iter + 1
    Sys.sleep(0.01)
    flush.console()
  }
  cat(paste0('[INFO] diff_beta: ', round(diff_beta, 7), "\n"))
  cat(paste0('[INFO] diff_theta: ', round(diff_theta, 7), "\n"))
  return(list(bhat=bhat, thetahat=thetahat))
}

get_sample_cov <- function(X, h_omega, h_beta){
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
  num_blocks,
  block_size,
  zeropos_list,
  seed=1
){
  set.seed(seed)
  sig.blocks <- vector(mode = "list", length = num_blocks)
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
  sig.est <- as.matrix(do.call(bdiag, sig.blocks))
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

