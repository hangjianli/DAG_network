lassoIdentTheta <- function(X, vers, zeropos, max.iter, lambda, lambda2){
  # Initialization ----------------------------------------------------------
  n <- dim(X)[1]
  p <- dim(X)[2]
  
  rho.est <- rep(1,p)
  rho.est[1] <- sqrt(n)/norm(as.matrix(X[,1]),type = "f")
  iter <- 1
  likeli <- rep(0,max.iter)
  rho_old <- rep(0, p)
  rho_diff_ <- phi_diff_ <- 1
  phi_old <- matrix(0,p,p)
  phi_diff <- rho_diff <- rep(0, max.iter-1)
  
  while(iter <= max.iter &  (phi_diff_ > 1e-4 || rho_diff_ > 1e-4)){
    phi.est <- matrix(0,p,p)
    set.seed(vers)
    first_col <- lars(x = as.matrix(cbind(0,X[,1])),
                  y = X[,2]*rho.est[2],
                  normalize = F,
                  intercept = F,
                  use.Gram = F
    )
    phi.est[,2] <- c((coef(first_col,s = (lambda/2),mode="lambda"))[-1], rep(0,p-1))
    #update phi--------------------------------------
    phi.temp <- foreach(i=2:(p-1),.combine = "cbind",.packages = c("glmnet")) %dopar% {
      lars.temp <- glmnet(x = X[,1:i],
                          y = X[,i+1]*rho.est[i+1],
                          alpha = 1,
                          lambda = lambda/(2*n),
                          intercept = F,
                          standardize = F, 
                          family = "gaussian",
                          thresh = 1e-10)
      c(as.numeric(lars.temp$beta), rep(0,p-i))
    }
    phi.est[,3:p] <- phi.temp
    #update rho------------------------------------
    rho.temp <- foreach(i = 2:p, .combine = "c") %dopar% {
      b_ <- -t(X[,i])%*%X%*%phi.est[,i]
      a_ <- t(X[,i])%*%X[,i]
      c_ <- -n
      (-b_ + sqrt(b_^2 - 4*a_*c_))/(2*a_)
    }
    rho.est[2:p] <-rho.temp
    #compute negative log likelihood---------------
    likeli.temp <- foreach(i=1:p, .combine = "c") %dopar% {
      -2*n*log(rho.est[i]) +
        norm(as.matrix(rho.est[i]*X[,i] - X%*%phi.est[,i]), type = "f")^2 +
        lambda*norm(as.matrix(phi.est[,i]),type = "1")
    }
    likeli[iter] <- sum(likeli.temp)
    #update iteration------------------------------
    rho_diff_ <- norm(as.matrix(rho_old - rho.est), type = "f") / sqrt(p)
    rho_old <- rho.est
    phi_diff_ <- norm((phi_old - phi.est), type = "f") / sqrt(p*(p-1)/2)
    phi_old <- phi.est
    phi_diff[iter] <- phi_diff_
    rho_diff[iter] <- rho_diff_
    iter <- iter + 1
    if(iter %% 10 == 0) cat("Iter: ", iter, "\n")
  }
  likeli <- round(likeli, 8)
  if(!all(diff(likeli[likeli!=0]) <= 0)){
    cat(likeli,"\n")
    stop("The likelihood for benchmark is not descreasing !!!")
  }
  
  b_est <- sweep(phi.est, 2, rho.est,'/')
  dimnames(b_est) <- list(dimnames(X)[[2]], dimnames(X)[[2]])
  s0_est <- sum(abs(b_est) > 1e-4)
  # saveRDS(b_est, paste0(n,p, "_", which(lambda.path == lambda[2]), "_b_bench.rds"))
  
  test.temp <- X%*%(diag(rho.est) - phi.est)
  test.res <- apply(test.temp,2,tcrossprod)
  S <- matrix(rowSums(test.res),n,n)
  S.scale <- S/p
  set.seed(vers+2)
  test <- glasso(s = S.scale,
                 thr = 1.0e-4,
                 rho = lambda2/p,
                 zero = zeropos,
                 penalize.diagonal = T)
  sig.est <- cov2cor(test$w)
  theta_est <- round(solve(sig.est),7)
  Theta_l0_lasso <- sum(abs(theta_est) > 1e-6)
  
  
  
  return(list(b_est = b_est,
              s0_est = s0_est,
              rho.est = rho.est,
              sig.est = sig.est,
              theta.est = theta_est,
              likeli_seq = likeli[likeli!=0],
              Theta_l0_lasso=Theta_l0_lasso))
}