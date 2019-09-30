dir.create(paste0("test",test_num))
setwd(paste0("test",test_num))

if(theta_type == "random"){
  sink(paste0("TwoStep-scaled-", test_num, "-p-",p,"-n-",n, "-theta_type-", theta_type,
              "-lambda2-",lambda2, "-magnitude-", b.magnitude,
              "-lambda1-", lambda[1],"-theta_sparsity-", sparse.prob, "-bProb-",
              prob_B, "-fix.zero-", fix.zero, '.txt'))
}else{
  sink(paste0("TwoStep-scaled-", test_num, "-p-",p,"-n-",n,
              # "-cliquesize-", clique_size,
              # "-epsilon-", epsilon,
              "-lambda2-",lambda2,
              # "-magnitude-", b.magnitude,
              "-lambda1-", lambda[1],
              # "-theta_sparsity-", sparse.prob, "-bProb-",
              # prob_B, "-fix.zero-", fix.zero,
              '.txt'))
}


# initialize --------------------------------------------------------------
omg.est <- rep(1,p)
theta_est3 <- diag(n)
choles <- chol(theta_est3) #upper triangular
# omg.est[1] <- sqrt(n/norm(as.matrix(choles%*%X[,1]),type = "f")^2)
beta.est <- matrix(0,p,p) 

# sacled lasso step --------------------------------------------------------------
set.seed(134)
# first <- scalreg(X = as.matrix(choles%*%cbind(X[,1])),
#                  y = as.matrix(choles%*%X[,2]*omg.est[2]),
#                  lam0 = lambda[2]/(2*n))
# beta.est[,2] <- c(first$coefficients, rep(0, p - 1))
# omg.est[2] <- first$hsigma
# 
# # update phi--------------------------------------
# phi.temp <- foreach(i=2:(p-1), .combine = "cbind", .packages = 'scalreg') %dopar% {
#   tst <- scalreg(X = as.matrix(choles)%*%X[,1:i],
#                  y = as.matrix(choles)%*%X[,i+1]*omg.est[i+1],
#                  lam0 = lambda[i+1]/(2*n))
#   c(tst$coefficients, rep(0,(p-i)), tst$hsigma)
# }
# beta.est[,3:p] <- phi.temp[-(p+1),]
# omg.est[3:p] <- phi.temp[p+1,]

phi.temp <- foreach(i=1:(p-1), .combine = "cbind", .packages = 'scalreg') %dopar% {
  if(i %% pp == 1) {
    first <- scalreg(X = as.matrix(choles%*%cbind(X[,i])),
                     y = as.matrix(choles%*%X[,i+1]*omg.est[i+1]),
                     lam0 = lambda[i+1]/(2*n))
    c(rep(0, pp * (i%/%pp)), c(first$coefficients, rep(0, pp - 1)), rep(0, pp * (copies - (i%/%pp) - 1)), first$hsigma)
  }else if(i %% pp != 0){
    tst <- scalreg(X = as.matrix(choles)%*%X[,(i - i%%pp + 1):i],
                   y = as.matrix(choles)%*%X[,i+1],
                   # *rho.est[i+1],
                   lam0 = lambda[i+1]/(2*n))
    c(rep(0, pp * (i%/%pp)), tst$coefficients,  rep(0, pp - i%%pp), rep(0, pp * (copies - (i%/%pp) - 1)), tst$hsigma)
  }
}

beta.est[, -seq(1, by =pp,  length.out = copies)] <- phi.temp[-(p+1),]
omg.est[seq(1, by =pp,  length.out = copies)] <- apply(X[,seq(1, by =pp,  length.out = copies), drop = F], 2, function(x) sqrt(n/norm(as.matrix(choles%*%x),type = "f")^2))
omg.est[-seq(1, by =pp,  length.out = copies)] <- phi.temp[p+1,]

test.temp <- X%*%(diag(1/omg.est) - sweep(beta.est, MARGIN = 2,STATS = omg.est, FUN = "/"))
test.res <- apply(test.temp, 2, tcrossprod)
S <- matrix(rowSums(test.res), n, n)
S.scale <- S/p

if(fix.zero){
  set.seed(245)
  test <- glasso(s = S.scale,
                 thr = 1.0e-4,
                 rho = lambda2/p,
                 zero = zeropos,
                 penalize.diagonal = T)
  # theta_est3 <- (test$wi + t(test$wi))/2  # manually symmetrize the inexact theta_est3 from glasso function
  sig.est <- cov2cor(test$w)
  theta_est3 <- round(solve(sig.est),7)
}else{
  set.seed(245)
  test <- glasso(s = S.scale,
                 thr = 1.0e-4,
                 rho = lambda2/p)
  sig.est <- cov2cor(test$w)
  # theta_est3 <- (test$wi + t(test$wi))/2  # manually symmetrize the inexact theta_est3 from glasso function
  # sig.est <- cov2cor(solve(theta_est3))
  theta_est3 <- round(solve(sig.est),7)
}

# total likelihood-----------------------------------


Theta_l0 <- sum(abs(theta_est3) > 1e-6)
B_l0 <- sum(abs(beta.est) > 1e-6)
b_est_heur <- beta.est
Theta_error_heur <- norm((theta - theta_est3), "f")/(n^2-dim(zeropos)[1])
Omega_error_heur <- sum((omg.est - omg)^2)/p
B_error_heuristic <- norm((b_est_heur - b), type = "f")/s0

saveRDS(b_est_heur, paste0(n,p, "_", which(lambda.path == lambda[2]), "_b_twostep.rds"))

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

cat("Omega estimates: ", 1/omg.est, "\n\n")
cat("True omega: ", omg, "\n\n")

sink()

