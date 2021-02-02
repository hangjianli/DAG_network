NLLglasso <- function(p, theta, S, lambda = 0){
  # calculate negative log likelihood value and its gradient-----------------
  NLL <- -p*log(det(theta)) + sum(diag(S%*%theta)) + sum(abs(theta))*lambda
  grad <- -p*solve(theta) + S
  # + lambda * sign(theta)
  return(list(value = NLL, gradient = grad))
}

hardthreshold <- function(A, thresh){
  # hardthresholding --------------------------------------------------------
  A[abs(A) < thresh] <- 0
  return(A)
}

# NLLelasnet <- function(phi,rho,theta, n, p, lam1, lam2, lam.e){
#   test.temp <- X%*%(diag(rho) - phi)
#   test.res <- apply(test.temp,2,tcrossprod)
#   S <- matrix(rowSums(test.res),n,n)
#   return(-2*n*sum(log(rho)) - p*log(det(theta)) + sum(diag(S%*%theta)) +
#            lam1*sum(abs(phi)) + lam.e*norm(phi,"f")^2 + lam2*sum(abs(theta)))
# }

sim_X <- function(vers, p, n, omg.sq, sig, b){
  #' simulate X from network DAG given its parameters
  #' 
  #' \code{sim_X} returns X and E matrices generated from the network DAG
  #' 
  #' @param 
  #' 
  set.seed(vers)
  eps_mat <- matrix(0, n, p)
  eps_mat[,1] <- mvrnorm(1, mu = rep(0, args$n), Sigma = omg.sq[1]*sig)
  eps_mat[,2] <- mvrnorm(1, mu = rep(0, args$n), Sigma = omg.sq[2]*sig)

  X <- matrix(0, n, p)
  
  X[,1] <- eps_mat[,1]
  X[,2] <- X[,1]*b[1,2] + eps_mat[,2]
    
  for(i in 3:p) {
    eps_mat[, i] <- mvrnorm(1, mu = rep(0, n), Sigma = omg.sq[i]*sig)
    X[,i] <- rowSums(sweep(X[,1:i-1], MARGIN = 2, b[1:i-1,i], "*")) + eps_mat[,i]
    if (i %% 50 == 0)
      cat("Getting ", i, "th column of X. \n" )
  }
  
  # psi_inv <- (diag(p) - b) %*% diag(1/omg.sq) %*% t(diag(p) - b)
  # L <- chol(psi_inv)
  # X <- X%*%t(L)
  # X <- eps_mat%*%solve(diag(p) - b)
  # X <- sweep(X, 2, STATS = apply(X,2,sd), FUN = "/")
  dimnames(X) <- list(NULL, as.character(1:p))
  return(list(X = X, eps_mat = eps_mat))
}

# sim_X_old <- function(vers, p, args, omg.sq, sig, b){
#   # simulate X from parameters ----------------------------------------------
#   set.seed(vers)
#   eps_mat <- matrix(0, args$n, p)
#   for(i in 1:p) {
#     eps_mat[, i] <- mvrnorm(1, mu = rep(0, args$n), Sigma = omg.sq[i]*sig)
#     cat("Getting ", i, "th column of X. \n" )
#   }
#   X <- eps_mat%*%solve(diag(p) - b)
#   return(list(X = X, eps_mat = eps_mat))
# }



calc_bench_errors <- function(estimands){
  n <- dim(estimands$b)[1]
  p <- dim(estimands$b)[2]
# Compare theta against identity mat, b against zero mat --------------------------------------
  theta_bench_err <- norm((diag(n) - estimands$theta), type = "f")/(n^2-dim(estimands$zeropos)[1])
  b_bench_err <- norm(estimands$b,type = "f")/estimands$s0
  
  return(list(theta_bench_err = theta_bench_err, 
              b_bench_err = b_bench_err))
}


get_average <- function(summary){
  # calculate average summary statistics ------------------------------------
  summary_ <- summary$summary
  estimand <- summary$estimands
  
  final <- list()
  final$Method <- summary_[[1]]$Method
  final$Threshold <- summary_[[1]]$Threshold
  final$lambda.pos <- Reduce("+", lapply(summary_, "[[", 3)) / length(summary_)
  final$P <- Reduce("+", lapply(summary_, "[[", 4)) / length(summary_)
  final$FP <- Reduce("+", lapply(summary_, "[[", 5)) / length(summary_)
  final$TP <- Reduce("+", lapply(summary_, "[[", 6)) / length(summary_)
  final$FDR <- final$FP / final$P
  final$SHD <- final$FP + (estimand$s0 - final$TP)
  final$JI <- final$TP / (final$P + estimand$s0 - final$TP)
  final$TE <- Reduce("+", lapply(summary_, "[[", 10)) / length(summary_)
  final$BE <- Reduce("+", lapply(summary_, "[[", 11)) / length(summary_)
  final <- as.data.frame(final)

  sink("average.txt", type="output")
  print(xtable(final, digits=c(0,0, 4,2,1,1,1,3,1,3,6,6), latex.environment = center),
        hline.after = c(-1, 0, 5, 10, 15, 20), include.rownames = F)
  sink()
  return(final)
}


ldlt_decomp <- function(A){
  # perform LDLT decomposition ----------------------------------------------
  L <- t(sweep(x = chol(A), MARGIN = 1, STATS = diag(chol(A)), FUN = "/"))
  D <- diag(diag(chol(A))^2)
  return(list(L=L,D=D))
}

Permute_X <- function(X, pi=NULL, seed = 10){
  # Given X, return X permutated-----------------------------------------
  if(!is.null(pi)){
    P <- as.matrix(sparseMatrix(pi, seq_along(pi), x = 1))
    X.p <- X %*% P
    return(list(order = pi, pMat = P, Xperm = X.p))
  }
  
  p <- dim(X)[2]
  set.seed(seed)
  random_P <- sample(p, p, replace = F)
  P <- as.matrix(sparseMatrix(random_P, seq_along(random_P), x = 1))
  X.p <- X %*% P
  dimnames(X.p) <- list(NULL, as.character(random_P))
  return(list(order = random_P, pMat = P, Xperm = X.p))
}

calc_interm_psi <- function(P, trueB, omg.sq){
  p <- dim(trueB)[1]
  psi_true <- solve(t(diag(p) - trueB)) %*% diag(omg.sq) %*% solve(diag(p) - trueB)
  psi_P <- t(P) %*% psi_true %*% P
  ldlt <- ldlt_decomp(round(psi_P, 10))
  B_p <-  diag(p) - t(solve(ldlt$L))
  return(list(B_p=B_p, psi_true = psi_true, psi_P = psi_P))
}


count_FN <- function(Xtrue, Xest) sum(abs(Xest[abs(Xtrue) > 1e-4]) < 1e-4)


# MLE ---------------------------------------------------------------------

dag_mle_estimation <- function(Bhat, Lhat, X){
  thresh <- 0
  Xnew <- as.matrix(Lhat%*%X)
  n <- dim(Xnew)[1]
  p <- dim(Bhat)[1]
  Bnew <- matrix(0,p,p)
  omg_new_sq <- rep(0,p)
  omg_new_sq[1] <- sum(Xnew[,1]^2)/n
  for(j in 2:p){
    nonzero_id <- which(abs(Bhat[,j]) > 1e-4)
    if(length(nonzero_id) >= n)
      stop(paste0("Too many nonzero entries in ", j, "th column"))
    if(length(nonzero_id) != 0){
      mj <- lm(Xnew[,j]~Xnew[,nonzero_id]-1)
      thresh <- thresh + sum(summary(mj)$coefficients[,2])
      Bnew[nonzero_id,j] <- coef(mj)
      omg_new_sq[j] <- sum(mj$residuals^2)/n
    }else{
      omg_new_sq[j] <- sum(Xnew[,j]^2)/n
    }
  }
  if(is.null(dimnames(Bhat))){
    dimnames(Bnew) <- list(as.character(1:p), as.character(1:p))  
  }else{
    dimnames(Bnew) <- dimnames(Bhat)
  }
  thresh <- thresh/sum(Bnew!=0)
  return(list(Bmle = Bnew,
              omgmlesq = omg_new_sq,
              threshold = thresh))
}


BIC_dag <- function(X, block_idx, bmle, omgmle, theta){
  n <- dim(X)[1]
  p <- dim(X)[2]
  s0 <- sum(abs(bmle) > 1e-4)
  t0 <- as.integer((sum(abs(theta) > 1e-4) - n) / 2)
  LX <- chol(theta)%*%X
  test.temp <- (LX - LX%*%bmle)%*%diag(1/sqrt(omgmle))
  test.res <- apply(test.temp,2,tcrossprod)
  S <- matrix(rowSums(test.res),n,n)
  tracetrm <- sum(diag(S))
  if(!is.null(block_idx)){
    thetatrm <- -p*sum(sapply(block_idx, function(x){log(det(as.matrix(theta[x,x])))}))  
  }else{
    thetatrm <- - p*log(det(theta))
  }
  negloglikelihood <- n*sum(log(omgmle)) + thetatrm + tracetrm
  # fn <- log(max(n,p)) * (s0 + t0) 
  fn <- 1 * (s0 + t0) 
  BIC <- negloglikelihood + fn
  return(list(BIC = BIC,
              negloglikelihood = negloglikelihood,
              s0 = s0,
              penalty=fn,
              thetatrm=thetatrm
  ))
}

# GES ---------------------------------------------------------------------
get_adjmat_from_fges <- function(edgelist, p, varnames){
  # return adjmatrix of cpdag
  myadjmat <- matrix(0, p, p)
  dimnames(myadjmat) <- list(varnames, varnames)
  for (i in 1:length(edgelist)) {
    par_name <- word(edgelist[i], 1)
    chil_name <- word(edgelist[i], -1)
    par_ind <- which(varnames == par_name)
    chil_ind <- which(varnames == chil_name)
    myadjmat[par_ind, chil_ind] <- 1
    if (grepl("<->", edgelist[i])  || grepl("---", edgelist[i]) ) {
      myadjmat[chil_ind, par_ind] <- 1
    }
  }
  return(myadjmat)
}


get_adjmat_from_pc <- function(res_pc, p){
  skl_list <- res_pc$edges #Show the result's edges
  skl_pc <- matrix(0,p,p)
  var_name <- as.character(sort(as.numeric(res_pc$nodes)))
  colnames(skl_pc) <- var_name
  for (i in 1:length(skl_list)) {
    par_name <- word(skl_list[i], 1)
    chil_name <- word(skl_list[i], -1)
    par_ind <- which((var_name == par_name) == 1)
    chil_ind <- which((var_name == chil_name) == 1)
    skl_pc[par_ind, chil_ind] <- 1
    if (grepl("<->", skl_list[i])  || grepl("---", skl_list[i]) ) {
      skl_pc[chil_ind, par_ind] <- 1
    }
  }
  rownames(skl_pc) <- var_name
  return(skl_pc)
}


generate_dag_pcalg <- function(p, s = 2*p / (p*(p-1)/2),
                               seed = 123){
  set.seed(seed)
  g <- pcalg::randomDAG(p, s, lB = -1, uB = 1)
  varnames <- g@nodes
  s0 <- numEdges(g)
  cat("[INFO] Number of edges in the generated DAG: ", s0, "\n")
  # save adjmat of true dag
  adjmat_true <- 1*(abs(as(g,"matrix")) > 0)
  #convert true dag to cpdag
  trueCPDAG <- dag2cpdag(g)
  adjmat_trueCPDAG <- as(trueCPDAG,"matrix")
  # generate b from dag
  btemp <- gen.B.from.btrue(p = p, B_true = adjmat_true)  
  b <- btemp$b
  return(list(b = b,
              s0 = s0,
              g=g,
              varnames = varnames,
              adjmat_true = adjmat_true,
              trueCPDAG = trueCPDAG,
              adjmat_trueCPDAG = adjmat_trueCPDAG))
}

calculate_total_likelihood_reparam <- function(X, rho.est, phi.est, theta.est, lam1, lam2){
  n <- dim(X)[1]
  p <- dim(X)[2]
  
  result <- -n*sum(log(rho.est^2)) - p*log(det(theta.est)) + 
    lam1*sum(abs(phi.est)) + lam2*sum(abs(theta.est))
  L <- chol(theta.est)
  naive <- 0
  for(i in 1:p){
    naive <- naive + norm(rho.est[i]*L%*%X[,i] - L%*%X%*%phi.est[,i],"f")^2 
  }
  result_lasso <- -n*sum(log(rho.est^2)) + lam1*sum(abs(phi.est)) + naive
  result_glasso <- - p*log(det(theta.est)) + lam2*sum(abs(theta.est)) + naive
  result <- result + naive
  return(list(result = result,
              result_lasso = result_lasso,
              result_glasso=result_glasso))
}

# calculate SHD -----------------------------------------------------------

compute_SHD_dag <- function(adj1, adj_true, s0){
  myshd <- 0
  if(dim(adj1)[1] != dim(adj_true)[1]){
    stop("Dimension doesn't match!")
  }
  n <- dim(adj1)[1]
  FP <- FN <- pnum <- 0
  for(i in 1:n){
    for(j in i:n){
      if(adj1[i,j] == 1){
        pnum <- pnum + 1
      }
      
      if((adj1[i,j] == 0) && (adj_true[i,j] == 1)){
        FN <- FN + 1
        next
      } 
      if((adj1[i,j] == 1) && (adj_true[i,j] == 0)){
        FP = FP + 1
      } 
    }
  }
  myshd = as.integer(FN + FP)
  TP = as.integer(s0 - FN)
  
  
  return(list(pnum = as.integer(pnum),
              FN = as.integer(FN),
              TP = TP,
              FDR = FP / (FP + TP),
              JI = TP / (pnum + s0 - TP),
              myshd = myshd))
}



compute_SHD_detail <- function(adj1, adj_true, s0){    
  # also computed total number of estiamted edges (directed + undirected)
  miss <- wrong_dir <- inv_miss <- ud <- du <- 0
  if(dim(adj1)[1] != dim(adj_true)[1]){
    stop("Dimension doesn't match!")
  }
  n <- dim(adj1)[1]
  pnum <- 0
  for(i in as.character(1:n)){
    for(j in as.character(i:n)){
      if (adj1[i,j] || adj1[j,i]){
        pnum = pnum + 1
      }
      if((adj1[i,j]+adj1[j,i]) == 0 && (adj_true[i,j] + adj_true[j,i]) != 0){
        miss = miss + 1
      }else if((adj1[i,j]+adj1[j,i]) != 0 && (adj_true[i,j] + adj_true[j,i]) == 0){
        inv_miss = inv_miss + 1
      }else if((adj1[i,j]+adj1[j,i]) == 1 && (adj_true[i,j] + adj_true[j,i]) == 1 && (adj1[i,j]!=adj_true[i,j])){
        wrong_dir = wrong_dir + 1
      }else if((adj_true[i,j] + adj_true[j,i]) == 2 && (adj1[i,j] + adj1[j,i]) == 1){
        ud = ud + 1
      }else if((adj_true[i,j] + adj_true[j,i]) == 1 && (adj1[i,j] + adj1[j,i]) == 2){
        du = du + 1
      }
      
    } 
  }
  myshd = miss + wrong_dir + inv_miss + ud + du
  reversed = wrong_dir + ud + du
  TP = s0 - miss
  return(list(pnum = pnum,
              FN = miss,
              # wrong_dir = wrong_dir,
              # FP = inv_miss,
              TP = TP,
              FDR = inv_miss / (inv_miss + TP),
              JI = TP / (pnum + s0 - TP),
              # ud = ud,
              # du = du,
              reversed = reversed,
              myshd = myshd))
}


compute_SHD_2step <- function(cpdag1, cpdag2, dag1,dag2){
  myshd <- 0
  if(dim(adj1)[1] != dim(adj_true)[1]){
    stop("Dimension doesn't match!")
  }
  n <- dim(adj1)[1]
  shieldmatrix <- matrix(0,n,n)
  for(i in 1:n){
    for(j in 1:n){
      if(dag1[i,j] == dag2[i,j]) shieldmatrix[i,j] = 1
    }
  }
  for(i in 1:n){
    for(j in i:n){
      if((adj1[i,j] != adj_true[i,j]) && (shieldmatrix[i,j] != 1)){
        myshd = myshd + 1
        next
      } 
      if((i!=j) && (adj1[j,i] != adj_true[j,i]) && (shieldmatrix[i,j] != 1)){
        myshd = myshd + 1
      } 
    }
  }
  return(myshd)
}



get_shd_for_several_methods <- function(kmainbic = NULL, kmaincor = NULL, 
                                        kbenchbic = NULL, kbenchcor = NULL,
                                        kflipbic = NULL, kflipcor = NULL,
                                        today_date, tt, gesshd, pcshd,
                                        adjmat_trueCPDAG, thresh = 0.2, estimands, ordered=F){
  
  REbench_bic <- REbench_cor <- REflip_bic <- theta_REflip_bic <- NULL
  s0 <- estimands$s0
  output_ordered <- output_unordered <- vector(mode = "list")
  
  if(ordered){
    if(!is.null(kbenchbic)){
      bestresult_bench_bic <- readRDS(file = paste0(format(today_date, "%Y-%m-%d"),
                                                    "-nsim-", tt,
                                                    "-vers-", kbenchbic, "-lam-", kbenchbic,
                                                    "-benchResult",".rds"))
      # adjmat_benchCPDAG <- bnstruct::dag.to.cpdag(1*(abs(bestresult_bench_bic$b_est) > thresh))
      adjmat_bench <- 1*(abs(bestresult_bench_bic$b_est) > thresh)
      shdXbench <- compute_SHD_dag(adj1 = adjmat_bench, adj_true = adjmat_trueCPDAG, s0)
      output_ordered$shdXbench = unlist(shdXbench)
    }
      
    if(!is.null(kbenchcor)){
      bestresult_bench_cor <- readRDS(file = paste0(format(today_date, "%Y-%m-%d"),
                                                    "-nsim-", tt,
                                                    "-vers-", kbenchcor, "-lam-", 
                                                    kbenchcor,
                                                    "-benchResult",".rds"))
      # adjmat_benchCPDAG_cor <- bnstruct::dag.to.cpdag(1*(abs(bestresult_bench_cor$b_est) > thresh))
      adjmat_bench_cor <- 1*(abs(bestresult_bench_cor$b_est) > thresh)
      shdXbenchCor <- compute_SHD_dag(adjmat_bench_cor, adjmat_trueCPDAG, s0)
      output_ordered$shdXbenchCor = unlist(shdXbenchCor)
    }
    if(!is.null(kmainbic)){
      bestresult_main_bic <- readRDS(file = paste0(format(today_date, "%Y-%m-%d"),
                                                   "-nsim-", tt,
                                                   "-vers-", kmainbic, "-lam-", kmainbic,
                                                   "-mainResult",".rds"))
      
      adjmat_main <- 1*(abs(bestresult_main_bic$b_est) > thresh)
      shdXmain <- compute_SHD_dag(adjmat_main, adjmat_trueCPDAG, s0)
      output_ordered$shdXmain=unlist(shdXmain)
    }    
    if(!is.null(kmaincor)){
      bestresult_main_cor <- readRDS(file = paste0(format(today_date, "%Y-%m-%d"),
                                                   "-nsim-", tt,
                                                   "-vers-", kmaincor, "-lam-", 
                                                   kmaincor,
                                                   "-mainResult",".rds"))
      adjmat_main_cor <- 1*(abs(bestresult_main_cor$b_est) > thresh)
      shdXmainCor <- compute_SHD_dag(adjmat_main_cor, adjmat_trueCPDAG, s0)
      output_ordered$shdXmainCor = unlist(shdXmainCor)
    }
  }else{
    if(!is.null(kmainbic)){
      bestresult_main_bic <- readRDS(file = paste0(format(today_date, "%Y-%m-%d"),
                                                   "-nsim-", tt,
                                                   "-vers-", kmainbic, "-lam-", kmainbic,
                                                   "-mainResult",".rds"))
      adjmat_mainCPDAG <- bnstruct::dag.to.cpdag(1*(abs(bestresult_main_bic$b_est) > thresh))
      shdXmain <- compute_SHD_detail(adjmat_mainCPDAG, adjmat_trueCPDAG, s0)
      output_unordered$shdXmain=unlist(shdXmain)
    }
    
    if(!is.null(kmaincor)){
      bestresult_main_cor <- readRDS(file = paste0(format(today_date, "%Y-%m-%d"),
                                                   "-nsim-", tt,
                                                   "-vers-", kmaincor, "-lam-", 
                                                   kmaincor,
                                                   "-mainResult",".rds"))
      adjmat_mainCPDAG_cor <- bnstruct::dag.to.cpdag(1*(abs(bestresult_main_cor$b_est) > thresh))
      shdXmainCor <- compute_SHD_detail(adjmat_mainCPDAG_cor, adjmat_trueCPDAG, s0)
      output_unordered$shdXmainCor = unlist(shdXmainCor)
    }
  }
  
  
  # relative erros ----------------------------------------------------------
  # bestresult_main_bic_MLE <- readRDS(file = paste0(format(today_date, "%Y-%m-%d"),
  #                                              "-nsim-", tt,
  #                                              "-vers-", kmainbic, "-lam-", kmainbic,
  #                                              "-mainMLEResult",".rds"))
  # 
  # 
  # 
  # bestresult_main_cor_MLE <- readRDS(file = paste0(format(today_date, "%Y-%m-%d"),
  #                                              "-nsim-", tt,
  #                                              "-vers-", kmaincor, "-lam-", 
  #                                              kmaincor,
  #                                              "-mainMLEResult",".rds"))
  # 
  # REmain_bic <- sqrt(sum((bestresult_main_bic_MLE$Bmle - estimands$b)^2) / sum(estimands$b^2))
  # REmain_cor <- sqrt(sum((bestresult_main_cor_MLE$Bmle - estimands$b)^2) / sum(estimands$b^2))
  # # REflip_cor <- sqrt(sum((bestresult_flip_cor$b_est - estimands$b)^2) / sum(estimands$b^2))
  # theta_REmain_bic <- sqrt(sum((bestresult_main_bic$theta_est - estimands$theta)^2) / sum(estimands$theta^2))
  # theta_REmain_cor <- sqrt(sum((bestresult_main_cor$theta_est - estimands$theta)^2) / sum(estimands$theta^2))
  # # theta_REflip_cor <- sqrt(sum((bestresult_flip_cor$theta_est - estimands$theta)^2) / sum(estimands$theta^2))
  # 
  # B_re <- list(REbench_bic=REbench_bic,
  #              REbench_cor=REbench_cor,
  #              REmain_bic=REmain_bic,
  #              REmain_cor=REmain_cor
  #              # REflip_bic=REflip_bic
  #              )
  # 
  # theta_re <- list(theta_REmain_bic=theta_REmain_bic,
  #                  theta_REmain_cor=theta_REmain_cor
  #                  # ,
  #                  # REflip_bic=theta_REflip_bic
  #                  )
  # saveRDS(B_re, file = "Bre.rds")
  # saveRDS(theta_re, file = "thetaRE.rds")
  # 

# output ------------------------------------------------------------------
  if (ordered){
    return(output_ordered)  
  }else{
    output_unordered$shdGES = gesshd
    output_unordered$shdPC = pcshd
    return(output_unordered)
  }
}


get_lam_path <- function(p, XX, rho.est, lambda.len, div=100){
  lambda.max <- rep(0, p)
  for(i in 1:(p-1)) lambda.max[i+1] <- norm(2*t(XX[,1:i])%*%(XX[,i+1]*rho.est[i+1]), type = "i")
  lambda.max <- max(lambda.max)
  lambda.path <- lseq(lambda.max/div, lambda.max/10, lambda.len)
  return(lambda.path)
}


# get_zeros ---------------------------------------------------------------
get_zeros <- function(theta, args, estimands){
  block_size = args$block_size
  num_blocks = ceiling(args$n/block_size)
  zeropos_temp <- vector(mode = "list")
  for(i in 1:num_blocks){
    cur_blocks <- theta[(1 + (i-1)*block_size): min(i*block_size, args$n), (1 + (i-1)*block_size): min(i*block_size, args$n)]
    zeropos_temp[[i]] <- which(abs(cur_blocks) < 1e-3, arr.ind = T)
  }
  return(zeropos_temp)
}

