# sim X -------------------------------------------------------------------
X_ <- sim_X(vers = 1, p = estimands$realp,
            args = args, omg.sq = estimands$omg.sq,
            sig = estimands$sig, b = estimands$b)
X <- X_$X
# generate lambda values --------------------------------------------------
XX <- X
lambda.path <- get_lam_path(estimands$realp, XX, rep(1,estimands$realp), 10, 100)
lambda <- lambda.path[5]
n <- dim(X)[1]
p <- dim(X)[2]
# some other settings -----------------------------------------------------
num_blocks <- ceiling(args$n/args$block_size)
block_size <- args$block_size
zeropos_list <- estimands$zeropos_list
