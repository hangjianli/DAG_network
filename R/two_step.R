dir.create(paste0("test",test_num))
setwd(paste0("test",test_num))

if(theta_type == "random"){
  sink(paste0("TwoStep-", test_num, "-p-",p,"-n-",n, "-theta_type-", theta_type,
              "-lambda2-",lambda2, "-magnitude-", b.magnitude,
              "-lambda1-", lambda[1],"-theta_sparsity-", sparse.prob, "-bProb-",
              prob_B, "-fix.zero-", fix.zero, '.txt'))
}else{
  sink(paste0("TwoStep-", test_num, "-p-",p,"-n-",n, "-cliquesize-", clique_size,
              "-epsilon-", epsilon, "-lambda2-",lambda2, "-magnitude-", b.magnitude,
              "-lambda1-", lambda[1],"-theta_sparsity-", sparse.prob, "-bProb-",
              prob_B, "-fix.zero-", fix.zero, '.txt'))
}


# initialize --------------------------------------------------------------
rho.est3 <- rep(1,p)
theta_est3 <- diag(n)
choles <- chol(theta_est3) #upper triangular
rho.est3[1] <- sqrt(n/norm(as.matrix(choles%*%X[,1]),type = "f")^2)
phi.est3 <- matrix(0,p,p) 


# lasso step --------------------------------------------------------------
set.seed(134)
first <- lars(x = as.matrix(choles%*%cbind(0,X[,1])), 
              y = as.matrix(choles%*%X[,2]*rho.est3[2]),
              normalize = F,
              intercept = F,
              use.Gram = F)
phi.est3[,2] <- c((coef(first, s = lambda[2]/2, mode = "lambda"))[-1], rep(0, p - 1))
# update phi--------------------------------------
phi.temp <- foreach(i=2:(p-1),.combine = "cbind",.packages = "lars") %dopar% {
  lars.temp <- lars(x = as.matrix(choles)%*%X[,1:i], 
                    y = as.matrix(choles)%*%X[,i+1]*rho.est3[i+1], 
                    normalize = F,
                    intercept = F,
                    use.Gram = F)
  c(coef(lars.temp,s = lambda[i+1]/2, mode = "lambda"), rep(0,(p-i)))
}
phi.est3[,3:p] <- phi.temp

#new method------------------------------------
rho.temp <- foreach(i = 2:p,.combine = "c") %dopar% {
  supp <- which(phi.est3[,i] != 0)
  if(length(supp) == 0) {
    e_j <- X[,i,drop = FALSE]
    norm(matrix(e_j),type = "f")/sqrt(n)
  }
  else {
    X_temp <- X[, supp, drop = FALSE]
    e_j <- (diag(n) - X_temp%*%solve(t(X_temp)%*%X_temp)%*%t(X_temp))%*%X[,i,drop = FALSE]
    norm(matrix(e_j),type = "f")/sqrt(n - dim(X_temp)[2])
  }
  
}
rho.est3[2:p] <- 1/rho.temp

# compute negative log likelihood------------
likeli.temp0 <- foreach(i = 1 : p, .combine = "c") %dopar% {
  -2*n*log(rho.est3[i]) +
    norm(as.matrix(rho.est3[i]*(choles%*%X[, i]) - choles%*%X%*%phi.est3[, i]), type = "f")^2 +
    lambda[i]*norm(as.matrix(phi.est3[,i]), "1")
}
NLL1 <- sum(likeli.temp0) - p*log(det(theta_est3)) + sum(abs(theta_est3))*lambda2
cat("NLL before glasso: ", NLL1, "\n\n")

# compute sample covariance-------------------------
test.temp <- X%*%(diag(rho.est3) - phi.est3)
test.res <- apply(test.temp, 2, tcrossprod)
S <- matrix(rowSums(test.res), n, n)
S.scale <- S/p

if(fix.zero){
  set.seed(245)
  test <- glasso(s = S.scale,
                 thr = 1.0e-4,
                 rho = lambda2/p,
                 zero = zeropos)
  theta_est3 <- (test$wi + t(test$wi))/2  # manually symmetrize the inexact theta_est3 from glasso function
  # sig.est <- cov2cor(solve(theta_est3))
  # theta_est3 <- round(solve(sig.est),7)
}else{
  set.seed(245)
  test <- glasso(s = S.scale,
                 thr = 1.0e-4,
                 rho = lambda2/p)
  theta_est3 <- (test$wi + t(test$wi))/2  # manually symmetrize the inexact theta_est3 from glasso function
  sig.est <- cov2cor(solve(theta_est3))
  theta_est3 <- round(solve(sig.est),7)
}

# total likelihood-----------------------------------
NLL2 <- -p*log(det(theta_est3)) + 
  sum(diag(S%*%theta_est3)) + 
  sum(abs(theta_est3))*lambda2 + 
  sum(sweep(abs(phi.est3),2, lambda, "*")) - 
  2*n*sum(log(rho.est3))

cat("NLL after glasso: ",  NLL2, "\n\n")

Theta_l0 <- sum(abs(theta_est3) > 1e-6)
B_l0 <- sum(abs(sweep(phi.est3, 2, rho.est3,'/')) > 1e-6)
 
Theta_error_heur <- norm((theta - theta_est3), "f")/(n^2-dim(zeropos)[1])
Omega_error_heur <- sum((1/rho.est3 - omg)^2)/p
b_est_heur <- sweep(phi.est3,2,rho.est3,'/')
saveRDS(b_est_heur, paste0(n,p, "_", which(lambda.path == lambda[2]), "_b_twostep.rds"))

B_error_heuristic <- norm((b_est_heur - b), type = "f")/(p*(p-1)/2)

num_FP_two <- sapply(threshold, function(x) sum(abs(b_est_heur[upper.tri(b_est_heur)][abs(b[upper.tri(b)]) < 1e-7]) > x))
num_TP_two <- sapply(threshold, function(x) sum(abs(b_est_heur[upper.tri(b_est_heur)][abs(b[upper.tri(b)]) > 1e-7]) > x))


hamming.dist3 <- sapply(threshold, function(x) sum(abs(b_est_heur[upper.tri(b_est_heur)][abs(b[upper.tri(b)]) > 1e-7]) <= x)) + num_FP_two



cat("Theta L0 norm: ", Theta_l0, "\n\n")
cat("Real theta l0: ", sum(abs(theta) > 1e-6), "\n\n")
cat("Beta L0 norm: ", B_l0, "\n\n")
cat("Real beta l0: ", sum(abs(b)>1e-6),"\n\n")

cat("Theta error: ", Theta_error_heur, "\n\n")
cat("benchmark theta: ", theta_bench, "  (compare to identity) " ,"\n\n") 

cat("Beta error: ", B_error_heuristic, "\n\n")
cat("benchmark beta: ", norm(b,type = "f")/(p*(p-1)/2), "\n\n")
cat("Omega error: ", Omega_error_heur, "\n\n")

cat("Omega estimates: ", 1/rho.est3, "\n\n")
cat("True omega: ", omg, "\n\n")

sink()

