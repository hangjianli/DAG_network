X_ <- sim_X(vers = 1, p = estimands$realp,
            args = args, omg.sq = estimands$omg.sq,
            sig = estimands$sig, b = estimands$b)
X <- X_$X
lambda2 <- 0.01
L <- diag(n)
L <- chol(estimands$theta)
image(as(cov(t(as.matrix(L%*%X))), "dgCMatrix"))
image(as(estimands$sig, 'dgCMatrix'))

XX <- X
lambda.path <- get_lam_path(estimands$realp, XX, rep(1,estimands$realp), 10, 100)
lambda <- lambda.path[5]
n <- dim(X)[1]
p <- dim(X)[2]

num_blocks <- ceiling(args$n/args$block_size)
block_size <- args$block_size
zeropos_list <- estimands$zeropos_list

# 
# newX <- sparsebnData(as.matrix(LX), type = "continuous")
# spasebn_noniid <- estimate.dag(newX)
# solution2_idx <- select.parameter(spasebn_noniid, newX)
# solution2 <- select(spasebn_noniid, index = solution2_idx)
# sparsebn_origin <- get.adjacency.matrix(solution2)
# shd2 <- compute_SHD_detail(sparsebn_origin, adjmat_trueCPDAG, estimands$s0)

iter <- 1
set.seed(11)
# A <- matrix(runif(n^2)*2-1, ncol=n) 
# estimands$sig <- cov2cor(t(A) %*% A)
# estimands$theta <- solve(estimands$sig)
# theta_est <- cov2cor(t(A) %*% A)
theta_est <- estimands$theta

max.iter <- 5000
total.likeli <- rep(0, max.iter)
theta_old <- matrix(0, n, n)
phi_old <- matrix(0,p,p)
theta_diff_ <- phi_diff_ <- 1

# developer reference --------------------------
theta_diff <- phi_diff <- rep(0, max.iter - 1)
estimands$theta[1:10,1:10]

while((iter <= max.iter) &&  (theta_diff_ > 1e-6 || phi_diff_ > 1e-6)){
  choles <- chol(theta_est)
  phi.est <- matrix(0,p,p) 
  first <- glmnet(x = as.matrix(choles%*%cbind(0,X[,1])),
                  y = as.matrix(choles%*%X[,2]),
                  alpha = 1,
                  intercept = F,
                  lambda = lambda/(2*n),
                  standardize = F,
                  family = "gaussian",
                  thresh = 1e-10)
  phi.est[,2] <- c(as.numeric(first$beta)[-1], rep(0, p - 1))
  #update phi--------------------------------------
  phi.temp <- foreach(i=2:(p-1),.combine = "cbind",.packages = "glmnet") %dopar% {
    lars.temp <- glmnet(x = as.matrix(choles)%*%X[,1:i],
                        y = as.matrix(choles)%*%X[,i+1],
                        alpha = 1,
                        intercept = F,
                        standardize = F,
                        lambda = lambda/(2*n),
                        family = "gaussian",
                        thresh = 1e-10)
    c(as.numeric(lars.temp$beta), rep(0, p - i))
  }
  phi.est[,3:p] <- phi.temp
  # compute sample covariance-------------------------
  S <- (X - X%*%phi.est)%*%t(X - X%*%phi.est)
  S.scale <- S/p
  # compute likelihood ------------------------------------------------------
  likeli2 <- calculate_total_likelihood_reparam(X = X, rho.est = rep(1,p), phi.est = phi.est,
                                                theta.est = theta_est, lam1 = lambda, 
                                                lam2 = lambda2)$result
  
  # if((iter > 1) && (likeli2 >= total.likeli[iter - 1])){
  #   cat("[INFO] Latest likelihood: ", likeli2, "\n")
  #   cat("[INFO] Previous likelihood: ", total.likeli[iter - 1], "\n")
  #   cat("[INFO] This is iteration ", iter, "\n")
  #   phi.est <- phi_old
  #   cat("[INFO] Lasso step is not decreasing. Algorithm converges. \n")
  #   warning(paste0("Lasso step is not decreasing at iter: ", iter, " Algorithm converges. \n"))
  #   break
  # }
  
  #estimate theta--------------------------------------
  sig.blocks <- vector(mode = "list", length = num_blocks)
  for(i in 1:num_blocks){
    zeros = estimands$zeropos_list[[i]]
    if(dim(zeros)[1] == 0)
      zeros = NULL
    temp_sig <- glasso(s = S.scale[(1+block_size*(i-1)) : min(block_size*i,n), 
                                   (1+block_size*(i-1)) : min(block_size*i,n)],
                       thr = 1.0e-15,
                       rho = lambda2/p,
                       zero = zeros,
                       penalize.diagonal = F)
    sig.blocks[[i]] <- temp_sig$w
  }
  sig.est <- as.matrix(do.call(bdiag, sig.blocks))
  theta_est <- solve(sig.est)
  
  # total likelihood-----------------------------------
  post_glasso_nll <- calculate_total_likelihood_reparam(X = X, rho.est = rep(1,p),
                                                        phi.est = phi.est, theta.est = theta_est,
                                                        lam1 = lambda, lam2 = lambda2)
    
  test <- -p*log(det(theta_est)) + sum(diag(S%*%theta_est)) +
    sum(abs(theta_est))*lambda2 + lambda*sum(abs(phi.est))
  
  assertthat::are_equal(test, post_glasso_nll$result)

  total.likeli[iter] <- post_glasso_nll$result
  
  phi_diff_ <- norm(as.matrix(phi_old - phi.est), "f") / sqrt(p*(p-1)/2)
  phi_old <- phi.est
  theta_diff_ <- norm((theta_old - theta_est), "f") / n
  theta_old <- theta_est
  phi_diff[iter] <- phi_diff_
  theta_diff[iter] <- theta_diff_
  if(iter %% 50 == 0){
    cat("Iter: ", iter, "\n",
        "nll: ", total.likeli[iter], "\n",
        "Diff (within):  ", total.likeli[iter] - likeli2, "\n",
        "Diff (btw): ", total.likeli[iter] - total.likeli[iter - 1] , "\n",
        "---------------------------------------------","\n")  
  }
  
  iter <- iter + 1

}


plot(total.likeli[total.likeli!=0], type = 'b', lwd = 1, pch = 16, col = 'red')

all(diff(total.likeli[total.likeli!=0]) < 0)

res1 <- phi.est
res12 <- theta_est
image(res12)
