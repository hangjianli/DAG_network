# Creating folder ---------------------------------------------------------
dir.create(paste0("test",test_num))
setwd(paste0("test",test_num))
if(theta_type == "random"){
  sink(paste0("Main-", test_num, "-p-",p,"-n-",n, "-theta_type-", theta_type,
              "-lambda2-",lambda2, "-magnitude-", b.magnitude,
              "-lambda1-", lambda[1],"-theta_sparsity-", sparse.prob, "-bProb-",
              prob_B, "-fix.zero-", fix.zero, '.txt'))
}else{
  sink(paste0("Main-", test_num, "-p-",p,"-n-",n,
              # "-epsilon-", epsilon,
              "-lambda2-",lambda2, 
              # "-magnitude-", b.magnitude,
              "-lambda1-", lambda[1],    
              # "-theta_sparsity-", sparse.prob, "-bProb-",
              # prob_B, "-fix.zero-", fix.zero,
              '.txt'))
}
###############################################################################

prod.omega <- first.term <- second.term <- third.term <- fourth.term <- fifth.term <- rep(0, max.iter)
Theta_l0 <- Theta_l1 <- B_l0 <- B_l1 <- Theta_error <- B_error <- Omega_error <- rep(0, max.iter)

iter2 <- 1
BLS_ind <- F
theta_est <- diag(n)
total.likeli <- rep(0, max.iter)
rho.est <- rep(1,p)
theta_old <- matrix(0, n, n)
theta_diff <- 1
rho_old <- rep(0,p)
rho_diff <- 1
phi_old <- matrix(0,p,p)
phi_diff <- 1
interp_ind <- 1
theta_difff_ <- phi_difff <- rho_difff <- rep(0, max.iter - 1)


while((iter2 <= max.iter) &&  (theta_diff > 1e-2 || rho_diff > 1e-2 || phi_diff > 1e-2)){
  choles <- chol(theta_est) #upper triangular
  rho.est[1] <- sqrt(n)/norm(as.matrix(choles%*%X[,1]),type = "f")
  phi.est <- matrix(0,p,p) 
  set.seed(37)
  
  phi.temp <- foreach(i = 1:(p-1), .combine = "cbind", .packages = c("glmnet", "lars")) %dopar% {
    if(i %% pp == 1) {
      first <- glmnet(x = as.matrix(choles %*% cbind(0,X[,i])),
                    y = as.matrix(choles %*% X[,i+1]*rho.est[i+1]),
                    alpha = 1,
                    lambda = lambda[i + 1]/(2*n),
                    intercept = F,
                    standardize = F)
      c(rep(0, pp * (i%/%pp)), as.numeric(first$beta)[-1], rep(0,pp-1), rep(0, pp * (copies - (i%/%pp) - 1)))
      
    }else if(i %% pp != 0){
      lars.temp2 <- glmnet(x = as.matrix(choles) %*% X[, (i - i%%pp + 1):i],
                           y = as.matrix(choles) %*% X[,i+1]*rho.est[i+1],
                           alpha = 1,
                           lambda = lambda[i+1]/(2*n),
                           intercept = F,
                           standardize = F)
      c(rep(0, pp * (i%/%pp)),
        as.numeric(lars.temp2$beta), rep(0, pp - i%%pp),
        rep(0, pp * (copies - (i%/%pp) - 1)))
    }
  }
  phi.est[, -seq(1, by =pp,  length.out = copies)] <- phi.temp
  # first <- glmnet(x = as.matrix(choles%*%cbind(0,X[,1])),
  #                 y = as.matrix(choles%*%X[,2]*rho.est[2]),
  #                 alpha = 1,
  #                 lambda = lambda[2]/(2*n),
  #                 intercept = F,
  #                 standardize = F)
  # phi.est[,2] <- c(as.numeric(first$beta)[-1], rep(0, p - 1))
  # #update phi--------------------------------------
  # phi.temp <- foreach(i=2:(p-1),.combine = "cbind",.packages = "glmnet") %dopar% {
  #   lars.temp <- glmnet(x = as.matrix(choles)%*%X[,1:i],
  #                       y = as.matrix(choles)%*%X[,i+1]*rho.est[i+1],
  #                       alpha = 1,
  #                       lambda = lambda[i+1]/(2*n),
  #                       intercept = F,
  #                       standardize = F)
  # 
  #   c(as.numeric(lars.temp$beta), rep(0,p-i))
  # }
  # phi.est[,3:p] <- phi.temp
  #update rho------------------------------------
  rho.temp <- foreach(i=2:p,.combine = "c") %dopar% {
    b_ <- -t(X[,i])%*%theta_est%*%X%*%phi.est[,i]
    a_ <- t(X[,i])%*%theta_est%*%X[,i]
    c_ <- -n
    (-b_ + sqrt(b_^2 - 4*a_*c_))/(2*a_)
  }
  rho.est[2:p] <- rho.temp
  #############################################
  ##### compute negative log likelihood ######
  #############################################
  # compute sample covariance-------------------------
  test.temp <- X%*%(diag(rho.est) - phi.est)
  test.res <- apply(test.temp,2,tcrossprod)
  S <- matrix(rowSums(test.res),n,n)
  S.scale <- S/p
 
  likeli2 <- -2*n*sum(log(rho.est)) + sum(diag(theta_est %*% S)) + 
    lambda[1]*sum(abs(phi.est)) - p*log(det(theta_est)) + sum(abs(theta_est))*lambda2
  
  if(iter2 > 0) cat("NLL Before lasso regression: ", total.likeli[iter2 - 1], "\n")
  cat("NLL after lasso regression: ", likeli2, "\n")
  if(iter2 > 1) cat("Lasso step is decreasing: ", likeli2 < total.likeli[iter2 - 1], "\n")
  if((iter2 > 1) && (likeli2 >= total.likeli[iter2 - 1])){
    warning("Lasso step is not decreasing!!")
    break
  }
  #estimate theta--------------------------------------
  if(fix.zero){
    set.seed(245)
    test <- glasso(s = S.scale,
                   thr = 1.0e-4,
                   rho = lambda2/p,
                   zero = zeropos,
                   penalize.diagonal = T)
    # theta_est <- (test$wi + t(test$wi))/2  # manually symmetrize the inexact theta_est from glasso function
    # sig.est <- cov2cor(solve(theta_est))
    sig.est <- cov2cor(test$w)
    theta_est <- round(solve(sig.est),7)
  }else{
    set.seed(245)
    test <- glasso(s = S.scale,
                   thr = 1.0e-4,
                   rho = lambda2/p)
    # theta_est <- (test$wi + t(test$wi))/2  # manually symmetrize the inexact theta_est from glasso function
    # sig.est <- cov2cor(solve(theta_est))
    sig.est <- cor2cor(test$w)
    theta_est <- round(solve(sig.est),5)
  }
  # total likelihood-----------------------------------
  post_glasso_nll <- -p*log(det(theta_est)) + sum(diag(S%*%theta_est)) +
    sum(abs(theta_est))*lambda2 + lambda[1]*sum(abs(phi.est)) - 2*n*sum(log(rho.est))
  
  missingpart <- sum(sweep(abs(phi.est),2, lambda, "*")) - 2*n*sum(log(rho.est))
  cat("Iter: ", iter2, "\n", "NLL after glasso but before BLS: ", post_glasso_nll, '\n')
  ########################################################################
  ########## backtracking line search   ##################################
  ########################################################################
  if ((likeli2 < post_glasso_nll) && (iter2 != 1))  {
    BLS_ind <- T
    cat("Line search activated! \n")
    for(interp_ind in 1:9){
      t <- 1/(2^interp_ind)
      theta_seq <- (1-t)*theta_old + t*theta_est
      sig_seq <- cov2cor(solve(theta_seq))
      theta_seq <- round(solve(sig_seq),7)
      eval_line <- NLLglasso(theta = theta_seq, p = p, S = S, lambda = lambda2)$value
      if(eval_line + missingpart < likeli2){
        cat("Line search stopped at t =", t, "\n")
        theta_est <- theta_seq
        total.likeli[iter2] <- eval_line + missingpart
        break
      }
      cat("t =", interp_ind, "NLL = ", eval_line + missingpart, "\n")
      if(interp_ind == 9){
        cat("Line search cannot improve the result. The algorithm converges. ")
        total.likeli[iter2] <- eval_line + missingpart
        break
      }
    }
  }else{
    total.likeli[iter2] <- post_glasso_nll
  }

  if(!isSymmetric(theta_est)) stop("Theta_est is not symmetric !!!") 
  
  cat("Is BLS activated ? ", BLS_ind, "\n")
  if(BLS_ind) cat("NLL after BLS: (should decrease)", "\n", total.likeli[iter2], "\n")
  cat("Decreasing after graphical lasso: ", likeli2 >= total.likeli[iter2], "\n",
      "Diff (within):  ", total.likeli[iter2] - likeli2, "\n",
      "Diff (btw): ", total.likeli[iter2] - total.likeli[iter2 - 1] , "\n",
      "===================================================","\n","\n")
  if(BLS_ind && interp_ind == 9){
    cat("Line search cannot improve theta. Algorithm terminates. ",'\n')
    theta_est <- theta_old
    break
  }
  
  phi_diff <- norm(as.matrix(phi_old - phi.est), "f")/sqrt(p*(p-1)/2)
  phi_old <- phi.est
  rho_diff <- norm(as.matrix(rho_old - rho.est), "f")/sqrt(p)
  rho_old <- rho.est
  theta_diff <- norm((theta_old - theta_est), "f")/n
  theta_old <- theta_est
  
  phi_difff[iter2] <- phi_diff
  rho_difff[iter2] <- rho_diff
  theta_difff_[iter2] <- theta_diff
  
  # compute some test statistic---------------------
  first.term[iter2] <- -2*n*sum(log(rho.est))
  second.term[iter2] <- -p*log(det(theta_est))
  third.term[iter2] <- sum(diag(S%*%theta_est))
  fourth.term[iter2] <- sum(sweep(abs(phi.est),2, lambda, "*"))
  prod.omega[iter2] <- prod(1/rho.est)
  fifth.term[iter2] <- lambda2*sum(abs(theta_est))
  
  ###########################
  ##### Count nonzeros ######
  ###########################
  
  Theta_l0[iter2] <- sum(abs(theta_est) > 1e-6)
  Theta_l1[iter2] <- sum(abs(theta_est))
  B_l0[iter2] <- sum(abs(sweep(phi.est, 2, rho.est,'/')) > 1e-6)
  B_l1[iter2] <- sum(abs(sweep(phi.est, 2, rho.est,'/')))
  
  ###########################
  ##### beta estimates ######
  ###########################
  
  b_est <- sweep(phi.est,2,rho.est,'/')
  
  ####################
  ##### errors #######
  ####################
  Omega_error[iter2] <- sum((1/rho.est - omg)^2)/p
  Theta_error[iter2] <- norm((theta - theta_est), type = "f")/(n^2-dim(zeropos)[1])
  B_error[iter2] <- norm(b_est - b, type = "f")/s0 
  iter2 <- iter2 + 1
  cat("Iter: ", iter2, "\n")
}


# Calculate summary statistics (P, TP, NP, SHD, JI, etc)--------------------------------------------
## A vector of length(threshold) of FP of b_est
num_FP_main <- sapply(threshold, function(x) sum(abs(b_est[upper.tri(b_est)][abs(b[upper.tri(b)]) < 1e-7]) > x))
## A vector of length(threshold) of TP of b_est
num_TP_main <- sapply(threshold, function(x) sum(abs(b_est[upper.tri(b_est)][abs(b[upper.tri(b)]) > 1e-7]) > x))
## A vector of length(threshold) of FP + FN for each threshold
hamming.dist <- sapply(threshold, function(x) sum(abs(b_est[upper.tri(b_est)][abs(b[upper.tri(b)]) > 1e-7]) <= x)) + num_FP_main

# save the b_est matrix
saveRDS(b_est, paste0(n,p, "_", which(lambda.path == lambda[2]), "_b_main.rds"))



cat("neg-sum-log-rho.est: ",first.term, "\n\n")
cat("neg-plog-det-theta_est: ", second.term, "\n\n")
cat("trace: ", third.term, "\n\n")
cat("lambda1 term: ", fourth.term, "\n\n")
cat("lambda2 term: ", fifth.term, "\n\n")

cat("Theta L0 norm: ", Theta_l0, "\n\n")
cat("Real theta l0: ", sum(abs(theta) > 1e-6), "\n\n")
cat("Beta L0 norm: ", B_l0, "\n\n")
cat("Real beta l0: ", s0 <- sum(abs(b)>1e-6),"\n\n")

cat("Theta error: ", Theta_error, "\n\n")
cat("True theta error: ", norm((theta - theta_est),type = "f")/(n^2-dim(zeropos)[1]), "\n")
cat("benchmark theta: ", theta_bench,"\n\n")

cat("Beta error: ", B_error, "\n\n")
cat("benchmark beta: ", b_bench, "\n\n")

cat("Omega error: ", Omega_error, "\n\n")
cat("Omega estimates: ", 1/rho.est, "\n\n")
cat("True omega: ", omg, "\n\n")

sink()

