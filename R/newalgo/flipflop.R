flipflop <- function(
  X,
  block_size,
  zeropos_list,
  block_idx=NULL,
  lambda2=100,
  lambda1=20,
  maxIter=10
){
  
  n <- dim(X)[1]
  p <- dim(X)[2]
  thetahat <- diag(n)
  num_blocks <- length(zeropos_list)
  
  psi_old <- matrix(0, p, p)
  theta_old <- matrix(0, n, n)
  theta_diff_ <- psi_diff_ <- 1
  iter <- 1
  loss_his <- vector(length = maxIter)
  if(n < p){
    theta_est <- solve(cov(t(X)))  
  }else{
    theta_est <- diag(n) 
  }
  
  total.likeli <- rep(0, maxIter)
  theta_diff <- psi_diff <- rep(0, maxIter - 1)
  
  while((iter <= maxIter) &&  (theta_diff_ > 1e-7 || psi_diff_>1e-7)){
    S.psi <- t(X)%*%theta_est%*%X / n
    set.seed(234)
    cat(paste0("===============================================================","\n"))
    cat("Start glasso for Psi\n")
    psi_est_inv <- estimate_theta(
      S = S.psi,
      p = n, 
      lambda2 = lambda1,
      block_size = p,
      zeropos_list = list(c())
    )
    #estimate theta--------------------------------------
    S.theta <- X%*%psi_est_inv%*%t(X) / p
    set.seed(1)
    cat("Start glasso for Theta\n")
    
    theta_est <- estimate_theta(
      S = S.theta, 
      p = p, 
      lambda2 = lambda2,
      block_size = block_size,
      zeropos_list = zeropos_list,
      block_idx = block_idx
    )
    
    # if (det(psi_est_inv) > 1e10 || det(psi_est_inv)< 1e-10){
    #   stop('The determinant of psi_est is unstable! \n')
    # }
    
    total.likeli[iter] <- -n*sum(log(eigen(psi_est_inv)$values))-
      p*sum(log(eigen(theta_est)$values)) + sum(diag(theta_est%*%X%*%psi_est_inv%*%t(X))) +
      lambda2*sum(abs(theta_est)) + lambda1*sum(abs(psi_est_inv))
    
    theta_diff_ <- norm((theta_old - theta_est), "f") / n
    psi_diff_ <- norm((psi_old - psi_est_inv), "f") / p
    theta_old <- theta_est
    psi_old <- psi_est_inv
    theta_diff[iter] <- theta_diff_
    psi_diff[iter] <- psi_diff_
    
    cat("[INFO] Iter: ", iter, "\n")
    cat("[INFO] theta_diff: ", theta_diff_, "psi_diff: ", psi_diff_, "\n")
    
    if(iter == 1){
      cat("[INFO] Saving estimates after one iteration. \n")
      psi_est <- solve(psi_est_inv)
      Ltemp <- t(chol(psi_est))
      L <- sweep(x = Ltemp, MARGIN = 1, STATS = diag(Ltemp), FUN = "/")
      bhat_1iter <- diag(p) - solve(t(L))
      dimnames(bhat_1iter) <- list(dimnames(X)[[2]], dimnames(X)[[2]])
      
      thetahat_1iter <- thetahat
    }
    iter <- iter + 1
    Sys.sleep(0.1)
    flush.console()
    
  }
  
  psi_est <- solve(psi_est_inv)
  Ltemp <- t(chol(psi_est))
  L <- sweep(x = Ltemp, MARGIN = 1, STATS = diag(Ltemp), FUN = "/")
  D <- diag(diag(Ltemp)^2)
  # norm(L%*%D%*%t(L) - psi_est,"f")
  bhat <- diag(p) - solve(t(L))
  dimnames(bhat) <- list(dimnames(X)[[2]], dimnames(X)[[2]])
  
  return(list(
    psi_est=psi_est,
    bhat=bhat, 
    omega2hat = as.numeric(diag(D)),
    bhat_1iter = bhat_1iter,
    thetahat=theta_est,
    thetahat_1iter=thetahat_1iter,
    losses=total.likeli[total.likeli!=0]
    )
  )
}
