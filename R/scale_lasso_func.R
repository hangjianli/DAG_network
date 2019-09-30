scale_lasso_init <- function(X, vers, lambda, lambda2, zeropos){
  
  n <- dim(X)[1]
  p <- dim(X)[2]
  
  omg.est <- rep(1,p)
  theta_est <- choles <- diag(n)
  omg.est[1] <- sqrt(n/norm(as.matrix(choles%*%X[,1]),type = "f")^2)
  beta.est <- matrix(0,p,p) 
  # sacled lasso step --------------------------------------------------------------
  set.seed(vers)
  first <- scalreg(X = as.matrix(choles%*%cbind(X[,1])),
                   y = as.matrix(choles%*%X[,2]*omg.est[2]),
                   lam0 = lambda/(2*n))
  beta.est[,2] <- c(first$coefficients, rep(0, p - 1))
  omg.est[2] <- first$hsigma
  # update phi--------------------------------------
  phi.temp <- foreach(i=2:(p-1), .combine = "cbind", .packages = 'scalreg') %dopar% {
    tst <- scalreg(X = as.matrix(choles)%*%X[,1:i],
                   y = as.matrix(choles)%*%X[,i+1]*omg.est[i+1],
                   lam0 = lambda/(2*n))
    c(tst$coefficients, rep(0,(p-i)), tst$hsigma)
  }
  beta.est[,3:p] <- phi.temp[-(p+1),]
  omg.est[3:p] <- phi.temp[p+1,]
  
  test.temp <- X%*%(diag(1/omg.est) - sweep(beta.est, MARGIN = 2, STATS = omg.est, FUN = "/"))
  test.res <- apply(test.temp, 2, tcrossprod)
  S <- matrix(rowSums(test.res), n, n)
  S.scale <- S/p
# Graphical lasso ---------------------------------------------------------
  set.seed(vers + 2)
  test <- glasso(s = S.scale,
                 thr = 1.0e-4,
                 rho = lambda2/p,
                 zero = zeropos,
                 penalize.diagonal = T)
  sig.est <- cov2cor(test$w)
  theta_est <- round(solve(sig.est),7)
  
  return(list(b_est = beta.est,
              omg.est = omg.est,
              sig.est = sig.est,
              theta.est = theta_est))
}
