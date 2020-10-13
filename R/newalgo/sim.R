source(file = 'R/loadpackages.R')
n <- 10
p <- 15

sqrt(log(p)/n)

omega <- sample(1:5, size = p, replace = T)

bchoice <- tolower(strsplit('pathfinder', NULL)[[1]][1])
omega

bname <- 'pathfinder'
btype <- 'discrete'
btrue <- try({constrB_from_BNrepo(name = bname, type = btype, ncopy = 1)})
pp = btrue$pp
b.temp <- gen.B.from.btrue(p = btrue$realp, seed = 1, btype = btype,
                           B_true = btrue$B_true, lower.thresh = 0.1, b.mag = 1)
b.temp$b %>% dim()

dim(permutations(52, 52))

blocks <- vector(mode = "list") # blocks of sigma
zeropos <- vector(mode = "list")
for(i in 1:num_blocks){
  if(i == num_blocks)
    block_size = n - (num_blocks-1)*block_size
  indices <- as.matrix(expand.grid(1:block_size, 1:block_size))
  blocks[[i]] <- diag(block_size)
  for(j in 1:nrow(indices)) blocks[[i]][indices[j,1], indices[j,2]] <- 0.7^(abs(diff(indices[j,]))/5)
  if(!is.positive.definite(blocks[[i]])) stop("Sigma is not positive definite !!!")
  theta <- round(solve(blocks[[i]]),5)
  zeropos[[i]] <- which(abs(as.matrix(theta)) < 1e-3, arr.ind = T)
}

sig <- do.call(bdiag, blocks)
theta<- round(solve(sig),5)
theta[abs(theta) < 1e-3] = 0
 
 
 
 

# new simulation ----------------------------------------------------------
 
X_ <- sim_X(vers = 1, 
            n = args$n,
            p = estimands$realp,
            omg.sq = estimands$omg.sq,
            sig = estimands$sig, b = estimands$b)
X <- X_$X
n <- dim(X)[1]
p <- dim(X)[2]

# omega -------------------------------------------------------------------
X_iid <- X[seq(1, args$n, by = args$block_size),]
X_iid %>% dim()
# estimate omega using i.i.d. samples
omega2_hat_iid <- estimate_omega(X_iid)
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


diag(estimands$theta)
diag(htheta)
