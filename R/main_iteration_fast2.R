main_iteration_fast2 <- function(
  X, vers, 
  estimands,
  args,
  max.iter = 50, 
  lambda = 0.01, 
  lambda2 = 0.1
){
  n <- dim(X)[1]
  p <- dim(X)[2]
  num_blocks <- ceiling(args$n/args$block_size)
  block_size <- args$block_size
  zeropos_list <- estimands$zeropos_list
  
  
  iter <- 1
  BLS_ind <- F
  interp_ind <- 1 # interpolation index
  
  theta_est <- diag(n)
  total.likeli <- rep(0, max.iter)
  rho.est <- rep(1,p)
  
  theta_old <- matrix(0, n, n)
  rho_old <- rep(0, p)
  phi_old <- matrix(0,p,p)
  theta_diff_ <- phi_diff_ <- rho_diff_ <- 1
  
  # developer reference --------------------------
  theta_diff <- phi_diff <- rho_diff <- rep(0, max.iter - 1)
  # ----------------------------------------------------------------
  
  while((iter <= max.iter) &&  (theta_diff_ > 1e-4 || rho_diff_ > 1e-4 || phi_diff_ > 1e-4)){
    choles <- chol(theta_est) #upper triangular
    rho.est[1] <- sqrt(n)/norm(as.matrix(choles%*%X[,1]),type = "f")
    phi.est <- matrix(0,p,p) 
    set.seed(vers)
    first <- glmnet(x = as.matrix(choles%*%cbind(0,X[,1])),
                    y = as.matrix(choles%*%X[,2]*rho.est[2]),
                    alpha = 1,
                    intercept = F,
                    lambda = lambda/(2*n),
                    standardize = F,
                    family = "gaussian",
                    thresh = 1e-10)
    # phi.est[,2] <- c(as.numeric(coef(first, s=lambda/(2*n)))[-(1:2)], rep(0, p - 1))
    phi.est[,2] <- c(as.numeric(first$beta)[-1], rep(0, p - 1))
    #update phi--------------------------------------
    phi.temp <- foreach(i=2:(p-1),.combine = "cbind",.packages = "glmnet") %dopar% {
      lars.temp <- glmnet(x = as.matrix(choles)%*%X[,1:i],
                          y = as.matrix(choles)%*%X[,i+1]*rho.est[i+1],
                          alpha = 1,
                          intercept = F,
                          standardize = F,
                          lambda = lambda/(2*n),
                          family = "gaussian",
                          thresh = 1e-10)
      # c(as.numeric(coef(lars.temp, s=lambda/(2*n)))[-1], rep(0,p-i))
      c(as.numeric(lars.temp$beta), rep(0, p - i))
    }
    phi.est[,3:p] <- phi.temp
    # Check number of nonzero is smaller than n
    if(!all(apply(phi.est,2,function(x) sum(x!=0)) <= n)){
      cat("[INFO] This is iteration ", iter, "\n")
      b_est <- sweep(phi.est,2,rho.est,'/')
      dimnames(b_est) <- list(dimnames(X)[[2]], dimnames(X)[[2]])
      s0_est <- sum(abs(b_est) > 1e-4)
      return(list(b_est = b_est, 
                  s0_est = s0_est,
                  theta_est = theta_est,
                  rho.est = rho.est, 
                  lik_seq = total.likeli[total.likeli!=0],
                  phi_diff = phi_diff,
                  rho_diff = rho_diff,
                  theta_diff = theta_diff))
      stop("Number of nonzero is larger than n! \n")
    }
    #update rho------------------------------------
    rho.temp <- foreach(i=2:p,.combine = "c") %dopar% {
      b_ <- -t(X[,i])%*%theta_est%*%X%*%phi.est[,i]
      a_ <- t(X[,i])%*%theta_est%*%X[,i]
      c_ <- -n
      (-b_ + sqrt(b_^2 - 4*a_*c_))/(2*a_)
    }
    rho.est[2:p] <- rho.temp
    # compute sample covariance-------------------------
    test.temp <- X%*%(diag(rho.est) - phi.est)
    test.res <- apply(test.temp,2,tcrossprod)
    S <- matrix(rowSums(test.res),n,n)
    S.scale <- S/p
    
    # compute likelihood ------------------------------------------------------
    likeli2 <- calculate_total_likelihood_reparam(X = X, rho.est = rho.est, phi.est = phi.est,
                                                  theta.est = theta_est, lam1 = lambda, 
                                                  lam2 = lambda2)$result
    
    if((iter > 1) && (round(likeli2, 10) >= round(total.likeli[iter - 1], 10))){
      cat("[INFO] Latest likelihood: ", likeli2, "\n")
      cat("[INFO] Previous likelihood: ", total.likeli[iter - 1], "\n")
      cat("[INFO] This is iteration ", iter, "\n")
      phi.est <- phi_old
      rho.est <- rho_old
      cat("[INFO] Lasso step is not decreasing. Algorithm converges. \n")
      warning(paste0("Lasso step is not decreasing at iter: ", iter, " Algorithm converges. \n"))
      break
    }
    
    #estimate theta--------------------------------------
    set.seed(vers + 2)
    sig.blocks <- vector(mode = "list", length = num_blocks)
    for(i in 1:num_blocks){
      zeros = estimands$zeropos_list[[i]]
      if(dim(zeros)[1] == 0)
        zeros = NULL
      temp_sig <- glasso(s = S.scale[(1+block_size*(i-1)) : min(block_size*i,n), (1+block_size*(i-1)) : min(block_size*i,n)],
                         thr = 1.0e-4,
                         rho = lambda2/p,
                         zero = zeros,
                         penalize.diagonal = T)
      sig.blocks[[i]] <- cov2cor(temp_sig$w) 
    }
    sig.est <- as.matrix(do.call(bdiag, sig.blocks))
    theta_est <- round(solve(sig.est), 5)
    
    
    # total likelihood-----------------------------------
    post_glasso_nll <- -p*log(det(theta_est)) + sum(diag(S%*%theta_est)) +
      sum(abs(theta_est))*lambda2 + lambda*sum(abs(phi.est)) - 2*n*sum(log(rho.est))
    missingpart <- sum(sweep(abs(phi.est),2, lambda, "*")) - 2*n*sum(log(rho.est))
    # backtracking line search----------------------------------
    if ((likeli2 < post_glasso_nll) && (iter != 1))  {
      BLS_ind <- T
      cat("Iter: ", iter, " Line search activated! \n")
      for(interp_ind in 1:9){
        t <- 1/(2^interp_ind)
        theta_seq <- (1-t)*theta_old + t*theta_est
        sig_seq <- cov2cor(solve(theta_seq))
        theta_seq <- round(solve(sig_seq),7)
        eval_line <- NLLglasso(theta = theta_seq, p = p, S = S, lambda = lambda2)$value
        if(eval_line + missingpart < likeli2){
          cat("Line search stopped at t =", t, "(0<t<1/2)", "\n")
          theta_est <- theta_seq
          total.likeli[iter] <- eval_line + missingpart
          break
        }
        cat("t =", interp_ind, ", NLL = ", eval_line + missingpart, "\n")
        if(interp_ind == 9){
          cat("Line search cannot improve the result. The algorithm converges. \n")
          total.likeli[iter] <- eval_line + missingpart
          break
        }
      }
    }else{
      total.likeli[iter] <- post_glasso_nll
    }
    
    if(!isSymmetric(round(theta_est), 7)) warning("Theta_est is not symmetric !!!") 
    if(BLS_ind) cat("NLL after BLS (should decrease): ", total.likeli[iter], "\n")
    if(BLS_ind && interp_ind == 9){
      cat("Algorithm terminates at iteration: ", iter, '\n')
      theta_est <- theta_old
      break
    }
    # if(round(likeli2, 10) < round(total.likeli[iter], 10)) stop("likelihood didn't increase after graphical Lasso! Something is wrong.")
    # if(likeli2 < total.likeli[iter]) 
      # stop("likelihood didn't increase after graphical Lasso! Something is wrong.")
    
    phi_diff_ <- norm(as.matrix(phi_old - phi.est), "f") / sqrt(p*(p-1)/2)
    phi_old <- phi.est
    rho_diff_ <- norm(as.matrix(rho_old - rho.est), "f") / sqrt(p)
    rho_old <- rho.est
    theta_diff_ <- norm((theta_old - theta_est), "f") / n
    theta_old <- theta_est
    phi_diff[iter] <- phi_diff_
    rho_diff[iter] <- rho_diff_
    theta_diff[iter] <- theta_diff_
    
    # compute some test statistic---------------------
    # Theta_l0[iter] <- sum(abs(theta_est) > 1e-6)
    # Theta_l1[iter] <- sum(abs(theta_est))
    # B_l0[iter] <- sum(abs(sweep(phi.est, 2, rho.est,'/')) > 1e-6)
    # B_l1[iter] <- sum(abs(sweep(phi.est, 2, rho.est,'/')))
    
    if(iter %% 10 == 0){
      cat("Iter: ", iter, "\n",
          "nll: ", total.likeli[iter], "\n",
          "Diff (within):  ", total.likeli[iter] - likeli2, "\n",
          "Diff (btw): ", total.likeli[iter] - total.likeli[iter - 1] , "\n",
          "---------------------------------------------","\n")
    }
    iter <- iter + 1
  }
  b_est <- sweep(phi.est,2,rho.est,'/')
  dimnames(b_est) <- list(dimnames(X)[[2]], dimnames(X)[[2]])
  s0_est <- sum(abs(b_est) > 1e-4)
  return(list(b_est = b_est, 
              s0_est = s0_est,
              theta_est = theta_est,
              rho.est = rho.est, 
              lik_seq = total.likeli[total.likeli!=0],
              phi_diff = phi_diff,
              rho_diff = rho_diff,
              theta_diff = theta_diff))
}