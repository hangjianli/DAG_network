run_bcd <- function(X, tol=1e-7){
  #' run BCD to estimate param
  #' 
  diff <- 1 
  n <- dim(X)[1]
  p <- dim(X)[2]
  thetahat <- diag(n)
  lambda.path <- get_lam_path(p, X, rep(1,p), 10, 100)
  while(diff > tol){
    bhat <- estimate_b(thetahat, X)
    S <- get_sample_cov(X, h_omega, bhat)
    thetahat <- estimate_theta(bhat, X)
    
  }
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
    # lambda <- cv.glmnet(as.matrix(choles)%*%XX[,1:i],
    #                     as.matrix(choles)%*%XX[,i+1])$lambda.min
    lars.temp <- glmnet(x = as.matrix(choles)%*%XX[,1:i],
                        y = as.matrix(choles)%*%XX[,i+1],
                        alpha = 1,
                        intercept = F,
                        lambda = lambda[i]/(2*n),
                        standardize = F,
                        family = "gaussian",
                        thresh = 1e-10)
    c(as.numeric(lars.temp$beta), rep(0, p - i))
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
                       thr = 1.0e-4,
                       rho = lambda2/p,
                       zero = zeros,
                       penalize.diagonal = T)
    sig.blocks[[i]] <- cov2cor(temp_sig$w) 
  }
  sig.est <- as.matrix(do.call(bdiag, sig.blocks))
  theta_est <- round(solve(sig.est), 5)
  return(theta_est)
}



estimate_omega <- function(X){
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

calculate_likelihood <- function(){
  
}
