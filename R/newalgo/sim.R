source(file = 'R/loadpackages.R')
n <- 10
p <- 15

# new simulation ----------------------------------------------------------


X_ <- sim_X(vers = 1, 
            n = args$n,
            p = estimands$realp,
            omg.sq = estimands$omg.sq,
            sig = estimands$sig, 
            b = estimands$b)
X <- X_$X
n <- dim(X)[1]
p <- dim(X)[2]

# omega -------------------------------------------------------------------
X_iid <- X[seq(1, args$n, by = args$block_size),]
X_iid %>% dim()
# estimate omega using i.i.d. samples
omega2_hat_iid <- estimate_omega_square(X_iid)
cbind(est=omega2_hat_iid, true=estimands$omg.sq)
plot(estimands$omg.sq, estimands$omg.sq - omega2_hat_iid)
plot(estimands$omg.sq, omega2_hat_iid)
norm(as.matrix(omega2_hat_iid - estimands$omg.sq), type = 'f') / estimands$realp

#estimate omega using entire data set
omega2_hat <- estimate_omega(X)
cbind(est=omega2_hat, true=estimands$omg.sq)
plot(estimands$omg.sq, estimands$omg.sq - omega2_hat)
plot(estimands$omg.sq, omega2_hat)
norm(as.matrix(omega2_hat - estimands$omg.sq), type = 'f') / estimands$realp

# beta estimate -----------------------------------------------------------

#lambda sequence

lambda.path <- get_lam_path(p, X, rep(1,p), 10, 100)
lambda.path

for(i in 1:10){
  test <- estimate_b(
    n = args$n,
    p = estimands$realp,
    X = X,
    theta_hat = estimands$theta,
    lambda= rep(lambda.path[i], estimands$realp)
  )
  plot(image(as(test, class(estimands$theta))))
}
# 2*seq(1:estimands$realp)

image(as(test, class(estimands$theta)))

image(as(estimands$b, class(estimands$theta)))

norm(estimands$b - diag(estimands$realp), type = 'f') / estimands$realp^2
norm(estimands$b - test, type = 'f') / estimands$realp^2

image(as(abs(estimands$b - test) > 100, class(estimands$theta))) 

test[test!=0] %>% summary()
estimands$b[estimands$b!=0] %>% summary()

(estimands$b * test >= 0 ) %>% sum()


hbeta <- estimate_b(
  n = args$n,
  p = estimands$realp,
  X = X,
  theta_hat = estimands$theta,
  lambda= rep(lambda.path[8], estimands$realp)
)
# theta estimate ----------------------------------------------------------

get_sample_cov <- function(X, h_omega, h_beta){
  rho_root <- 1 / h_omega
  test.temp <- X%*%(diag(rho_root) - sweep(h_beta, MARGIN = 2, FUN = '*', STATS = rho_root))
  test.res <- apply(test.temp,2,tcrossprod)
  S <- matrix(rowSums(test.res), n, n)
  S.scale <- S/p
  return(S.scale)
}

S <- get_sample_cov(X, omega2_hat_iid, hbeta)

num_blocks <- ceiling(args$n/args$block_size)
block_size <- args$block_size
zeropos_list <- estimands$zeropos_list

htheta <- estimate_theta(S = S, 
                         p =p,
                         lambda2 = 1,
                         num_blocks = num_blocks,
                         block_size = block_size,
                         zeropos_list = zeropos_list)

image(estimands$theta)
image(as(htheta, class(estimands$theta)))
estimands$theta[1:10,1:10]
diag(estimands$theta)
diag(htheta)


# test bcd ----------------------------------------------------------------
lambda.path <- get_lam_path(p, X, rho.est = rep(1,p), 10, 100)
lambda.path
BICscores_main <- BICscores_flipflop <- BICscores_bench <- rep(0, 10)
minrowcor_main <- minrowcor_bench <- minrowcor_flip <- rep(0, 10)

k <- 7
for(k in 1:length(lambda.path)){
  res <- run_bcd(
    X = X, 
    block_size = args$block_size,
    zeropos_list = estimands$zeropos_list,
    lambda1 = lambda.path[k],
    lambda2 = .001,
    estimands = estimands, 
    tol = 1e-7)
  
  cor_est <- cor(t(chol(res$thetahat)%*%(X - X%*%res$bhat)))
  minrowcor_main[k] <- sum(abs(cor_est[upper.tri(cor_est)]))
  mle_result <- dag_mle_estimation(X = X, Bhat = res$bhat, Lhat = chol(res$thetahat))
  # saveRDS(mle_result, 
  #         file = paste0(format(today_date, "%Y-%m-%d"), 
  #                       "-nsim-", sim,
  #                       "-vers-", k, "-lam-", 
  #                       k,"-mainMLEResult",".rds")
  #         )
  BIC_result <- BIC_dag(
    X = X,
    bmle = mle_result$Bmle, 
    omgmle = mle_result$omgmlesq,
    theta = res$thetahat
  )
  BICscores_main[k] <- BIC_result$BIC
  # saveRDS(BICscores_main, "BICscores_main.rds")
  # saveRDS(minrowcor_main, "minrowcor_main.rds")
  
}

bestk_bic_main <- which.min(BICscores_main[BICscores_main!=0])
bestk_cor_main <- which.min(minrowcor_main[minrowcor_main!=0])
# bestk_bic_bench <- which.min(BICscores_bench[BICscores_bench!=0])
# bestk_cor_bench <- which.min(minrowcor_bench[minrowcor_bench!=0])

plot(BICscores_main[BICscores_main!=0], pch = 16, type = "b")
plot(minrowcor_main[minrowcor_main!=0], pch = 16, type = "b")

image(as(mle_result$Bmle, class(estimands$theta)))
image(as(res$bhat, class(estimands$theta)))
image(as(res$thetahat, class(estimands$theta)))

bhat_adj <- 1*(abs(res$bhat) > 0.0001)
bstar_adj <- 1*(abs(estimands$b) > 0.0)

thetahat_adj <- 1*(abs(res$thetahat) > 0)
thetastar_adj <- 1*(abs(estimands$theta) > 0)

image(as(bhat_adj, class(estimands$theta)))
image(as(bstar_adj, class(estimands$theta)))

bhat_stats <- compute_SHD_dag(adj1 = bhat_adj, adj_true = bstar_adj, estimands$s0) %>% unlist()
bhat_stats

# baseline ----------------------------------------------------------------

for(k in 1:length(lambda.path)){
  res_baseline <- run_bcd(
    X = X, 
    baseline_flag = T,
    block_size = args$block_size,
    zeropos_list = estimands$zeropos_list,
    lambda1 = lambda.path[k],
    lambda2 = 1,
    estimands = estimands, 
    tol = 1e-7)
  
  cor_est <- cor(t(chol(res_baseline$thetahat)%*%(X - X%*%res_baseline$bhat)))
  minrowcor_bench[k] <- sum(abs(cor_est[upper.tri(cor_est)]))
  mle_result_bench <- dag_mle_estimation(X = X, Bhat = res_baseline$bhat, Lhat = chol(res_baseline$thetahat))
  # saveRDS(mle_result, 
  #         file = paste0(format(today_date, "%Y-%m-%d"), 
  #                       "-nsim-", sim,
  #                       "-vers-", k, "-lam-", 
  #                       k,"-mainMLEResult",".rds")
  #         )
  BIC_result <- BIC_dag(
    X = X,
    bmle = mle_result_bench$Bmle, 
    omgmle = mle_result_bench$omgmlesq,
    theta = res_baseline$thetahat
  )
  BICscores_bench[k] <- BIC_result$BIC
  # saveRDS(BICscores_main, "BICscores_main.rds")
  # saveRDS(minrowcor_main, "minrowcor_main.rds")
  
}

allShdS <- get_shd_for_several_methods(
  kmainbic = bestk_bic_main,
  kbenchbic = bestk_bic_bench,
  kmaincor = bestk_cor_main,
  kbenchcor = bestk_cor_bench,
  today_date = '',
  tt = sim,
  gesshd = NULL, 
  pcshd = NULL,
  adjmat_trueCPDAG = adjmat_trueCPDAG,
  thresh = 0.1,
  estimands = estimands,
  ordered = T
)

image(as(res_baseline$bhat, class(estimands$theta)))
image(as(res_baseline$thetahat, class(estimands$theta)))


bhat_adj_baseline <- 1*(abs(res_baseline$bhat) > 0.02)
image(as(bhat_adj_baseline, class(estimands$theta)))

bhat_stats_baseline <- compute_SHD_dag(adj1 = bhat_adj_baseline, adj_true = bstar_adj, estimands$s0) %>% unlist()
bhat_stats_baseline


compute_SHD_dag(adj1 = matrix(0, p, p), adj_true = bstar_adj, estimands$s0) %>% unlist()

