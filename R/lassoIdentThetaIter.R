dir.create(paste0("test",test_num))
setwd(paste0("test",test_num))

rho.est2 <- rep(1,p)
rho.est2[1] <- sqrt(n)/norm(as.matrix(X[,1]),type = "f")
iter <- 1
likeli <- rep(0, nrow = max.iter)
rho_old2 <- rep(0,p)
rho_diff2 <- 1
phi_old2 <- matrix(0,p,p)
phi_diff2 <- 1
phi_diff_ <- rho_diff_ <- rep(0, max.iter2 - 1)

while(iter < max.iter2 &  (phi_diff2 > 1e-4 || rho_diff2 > 1e-4)){
  phi.est2 <- matrix(0,p,p)
  set.seed(134)
  # first <- lars(x = as.matrix(cbind(0,X[,1])),
  #               y = X[,2]*rho.est2[2],
  #               normalize = F,
  #               intercept = F,
  #               use.Gram = F
  # )
  # phi.est2[,2] <- c((coef(first,s = (lambda[2]/2),mode="lambda"))[-1], rep(0,p-1))
  # #update phi--------------------------------------
  # phi.temp <- foreach(i=2:(p-1),.combine = "cbind",.packages = c("glmnet")) %dopar% {
  #   lars.temp2 <- glmnet(x = X[,1:i],
  #                       y = X[,i+1]*rho.est2[i+1],
  #                       alpha = 1,
  #                       lambda = lambda[i+1]/(2*n),
  #                       intercept = F,
  #                       standardize = F)
  # 
  # 
  #   # c(lars.temp$B0, rep(0,(p-i)))
  #   c(as.numeric(lars.temp2$beta), rep(0,p-i))
  # }
  # phi.est2[,3:p] <- phi.temp

  phi.temp <- foreach(i = 1:(p-1),.combine = "cbind", .packages = c("glmnet", "lars")) %dopar%{
    if(i %% pp == 1) {
      first <- lars(x = as.matrix(cbind(0,X[,i])),
                                  y = X[,i+1]*rho.est2[i+1],
                                  normalize = F,
                                  intercept = F,
                                  use.Gram = F
                    )
      c(rep(0, pp * (i%/%pp)), (coef(first,s = (lambda[i+1]/2),mode="lambda"))[-1], rep(0,pp-1), rep(0, pp * (copies - (i%/%pp) - 1)))
    }else if(i %% pp != 0){
      lars.temp2 <- glmnet(x = X[, (i - i%%pp + 1):i],
                          y = X[,i+1]*rho.est2[i+1],
                          alpha = 1,
                          lambda = lambda[i+1]/(2*n),
                          intercept = F,
                          standardize = F)
    c(rep(0, pp * (i%/%pp)),
      as.numeric(lars.temp2$beta), rep(0, pp - i%%pp),
      rep(0, pp * (copies - (i%/%pp) - 1)))
    }
  }
  phi.est2[, -seq(1, by =pp,  length.out = copies)] <- phi.temp
  #update rho------------------------------------
  rho.temp <- foreach(i = 2:p, .combine = "c") %dopar% {
    b_ <- -t(X[,i])%*%X%*%phi.est2[,i]
    a_ <- t(X[,i])%*%X[,i]
    c_ <- -n
    (-b_ + sqrt(b_^2 - 4*a_*c_))/(2*a_)
  }
  rho.est2[2:p] <-rho.temp
  #compute negative log likelihood---------------
  likeli.temp <- foreach(i=1:p, .combine = "c") %dopar% {
    -2*n*log(rho.est2[i]) +
      norm(as.matrix(rho.est2[i]*X[,i] - X%*%phi.est2[,i]), type = "f")^2 + 
      lambda[i]*norm(as.matrix(phi.est2[,i]),type = "1")
  }
  likeli[iter] <- sum(likeli.temp)
  #update iteration------------------------------
  rho_diff2 <- norm(as.matrix(rho_old2 - rho.est2), type = "f")
  rho_old2 <- rho.est2
  phi_diff2 <- norm((phi_old2 - phi.est2), type = "f")
  phi_old2 <- phi.est2
  phi_diff_[iter] <- phi_diff2
  rho_diff_[iter] <- rho_diff2
  iter <- iter + 1
  cat("Iter: ", iter, "\n")
}

b_est2 <- sweep(phi.est2, 2, rho.est2,'/')
saveRDS(b_est2, paste0(n,p, "_", which(lambda.path == lambda[2]), "_b_bench.rds"))
num_FP_bench <- sapply(threshold, function(x) sum(abs(b_est2[upper.tri(b_est2)][abs(b[upper.tri(b)]) < 1e-7]) > x))
num_TP_bench <- sapply(threshold, function(x) sum(abs(b_est2[upper.tri(b_est2)][abs(b[upper.tri(b)]) > 1e-7]) > x))
hamming.dist2 <- sapply(threshold, function(x) sum(abs(b_est2[upper.tri(b_est2)][abs(b[upper.tri(b)]) > 1e-7]) <= x)) + num_FP_bench # FN + FP
B_error_bench <- norm((b_est2 - b), type = "f")


test.temp <- X%*%(diag(rho.est2) - phi.est2)
test.res <- apply(test.temp,2,tcrossprod)
S <- matrix(rowSums(test.res),n,n)
S.scale <- S/p
set.seed(245)
test <- glasso(s = S.scale,
               thr = 1.0e-4,
               rho = lambda2/p,
               zero = zeropos,
               penalize.diagonal = T)
sig.est <- cov2cor(test$w)
theta_est <- round(solve(sig.est),7)
Theta_l0_lasso <- sum(abs(theta_est) > 1e-6)
Theta_error_larsvar <- norm((theta - theta_est), "f")/(n^2-dim(zeropos)[1])

