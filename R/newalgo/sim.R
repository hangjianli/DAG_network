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
omega2_hat <- estimate_omega_square(X)
cbind(est=omega2_hat, true=estimands$omg.sq)
plot(estimands$omg.sq, estimands$omg.sq - omega2_hat)
plot(estimands$omg.sq, omega2_hat)
norm(as.matrix(omega2_hat - estimands$omg.sq), type = 'f') / estimands$realp

# beta estimate -----------------------------------------------------------

#lambda sequence
lambda.path <- get_lam_path(p, X, rep(1,p), 10, 40)
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
  lambda= rep(lambda.path[14], estimands$realp)
)

bhat_adj_old <- 1*(abs(old_b$b_est) > 0.1)
bhat_adj <- 1*(abs(hbeta) > 0.1)
bstar_adj <- 1*(abs(estimands$b) > 0)

image(as(bhat_adj, class(estimands$theta)))
image(as(bhat_adj_old, class(estimands$theta)))
image(as(bstar_adj, class(estimands$theta)))

norm(hbeta - estimands$b, 'f')
norm(old_b$b_est - estimands$b, 'f')
bhat_stats <- compute_SHD_dag(adj1 = bhat_adj, adj_true = bstar_adj, estimands$s0) %>% unlist()
bhat_stats
# theta estimate ----------------------------------------------------------
omega2_hat_iid <- rep(1, estimands$realp)
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
norm(htheta - estimands$theta, 'f')^2 / args$n^2

# test bcd ----------------------------------------------------------------
rm(list = ls())
setwd("~/Documents/research/dag_network/")
source("R/loadpackages.R")
source("R/helper_funcs.R")
source("R/gen_params_funcs.R")
source("R/newalgo/functions.R")
# generate specs ----------------------------------------------------------
args <- args_for_parameter()  
dir.create(path = paste0('output/', args$setting))
saveRDS(args, file = paste0('output/', args$setting, "/args.rds"))
estimands <- generate_parameters(args = args, seed = 1)
image(estimands$theta)
image(as(estimands$b, class(estimands$theta)))
estimands$b[estimands$b!=0] %>% as.numeric() %>% hist(breaks=20)
saveRDS(estimands, file = paste0('output/', args$setting, "/estimands.rds"))
# run simulation ----------------------------------------------------------
setwd("~/Documents/research/dag_network")
sim_newalgo_ordered(args, estimands, start_sim=1, end_sim=args$num_sim, lamLen=15)
process_output_ordered(simID='001')
get_all_shd_ordered(simID, estimands, args$num_sim)






