flipflop <- function(X, n, p,
                     args,
                     zeropos_psi= NULL, 
                     max.iter = 15, 
                     lambda, lambda2,
                     zeropos_list){
  
  num_blocks <- ceiling(args$n/args$block_size)
  block_size <- args$block_size
  iter <- 1
  theta_est <- psi_est_inv <- diag(n)
  total.likeli <- rep(0, max.iter)
  theta_old <- matrix(0, n, n)
  psi_old <- matrix(0,p,p)
  theta_diff_ <- psi_diff_ <- 1
  theta_diff <- psi_diff <- rep(0, max.iter - 1)
  
  while((iter <= max.iter) &&  (theta_diff_ > 1e-4 || psi_diff_>1e-4)){
    S.psi <- t(X)%*%theta_est%*%X / n
    # S.psi <- t(X)%*%estimands$theta%*%X / n
    set.seed(234)
    cat("Start glasso for Psi\n")
    first <- glasso(s = as.matrix(S.psi),
                    thr = 1.0e-4,
                    rho = 0.001 / n,
                    zero = zeropos_psi[1:500,],
                    # wi.init = psi_est_inv,
                    penalize.diagonal = T,
                    maxit = 5, trace = T)
    # sum(psi_est_inv!=0)
    # image(as(psi_est_inv, "dgCMatrix"))
    # psi.inv <- first$w
    psi_est_inv <- round(first$wi, 7)
    #estimate theta--------------------------------------
    S.theta <- X%*%psi_est_inv%*%t(X) / p
    set.seed(245)
    
    cat("Start glasso for Theta\n")
    set.seed(1)
    sig.blocks <- vector(mode = "list", length = num_blocks)
    for(i in 1:num_blocks){
      zeros = zeropos_list[[i]]
      if(dim(zeros)[1] == 0)
        zeros = NULL
      temp_sig <- glasso(s = S.theta[(1+block_size*(i-1)):min(block_size*i,n), 
                                     (1+block_size*(i-1)):min(block_size*i,n)],
                         thr = 1.0e-3,
                         rho = lambda2/p,
                         zero = zeros,
                         penalize.diagonal = T,
                         maxit = 1)
      sig.blocks[[i]] <- cov2cor(temp_sig$w)  
    }
    sig.est <- as.matrix(do.call(bdiag, sig.blocks))
    theta_est <- round(solve(sig.est), 5)
    
    # test <- as(resultflip$b_est, "dgCMatrix")
    # image(test)
    
    # second <- glasso(s = S.theta,
    #                thr = 1.0e-5,
    #                rho = lambda2/p,
    #                zero = zeropos_theta,
    #                penalize.diagonal = T,
    #                maxit = 100)
    # sig.est <- cov2cor(second$w)
    # theta_est <- round(solve(sig.est),7)
    
    total.likeli[iter] <- -n*log(det(psi_est_inv)) -
      p*log(det(theta_est)) + sum(diag(theta_est%*%X%*%psi_est_inv%*%t(X))) +
      lambda2*sum(abs(theta_est)) + lambda*sum(abs(psi_est_inv))
    
    theta_diff_ <- norm((theta_old - theta_est), "f") / n
    psi_diff_ <- norm((psi_old - psi_est_inv), "f") / p
    theta_old <- theta_est
    psi_old <- psi_est_inv
    theta_diff[iter] <- theta_diff_
    psi_diff[iter] <- psi_diff_
    
    # if(iter %% 5 == 0){
    cat("[INFO] Iter: ", iter, "\n")
    cat("[INFO] theta_diff: ", theta_diff_, "psi_diff: ", psi_diff_, "\n")
    # } 
    iter <- iter + 1
  }
  
  psi_est <- solve(psi_est_inv)
  
  Ltemp <- t(chol(psi_est))
  L <- t(sweep(x = Ltemp, MARGIN = 1, STATS = diag(Ltemp), FUN = "/"))
  D <- diag(diag(Ltemp)^2)
  # norm(L%*%D%*%t(L) - psi_est,"f")

  b_est <- diag(p) - solve(t(L))
  dimnames(b_est) <- list(dimnames(X)[[2]], dimnames(X)[[2]])
  # heatmap.2(abs(b_est), dendrogram = "none", Rowv = F, Colv = F, trace = "none")
  
  return(list(psi_est = psi_est, 
              theta_est = theta_est,
              b_est = b_est,
              # rho.est = rho.est,
              likeli_seq=total.likeli,
              totaliter = iter,
              theta_diff = theta_diff_,
              psi_diff = psi_diff_))
}
