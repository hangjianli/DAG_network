#######################################
###### Generates omega, B, Theta ###### 
#######################################


# generate Omega-------------------
set.seed(23)
omg <- runif(p) # omega's between 0,1
omg[1] <- 1
omg.sq <- omg^2
# generate Beta---------------------
set.seed(482) # generate B, taking values in [-thresh, +thresh]
# b <- matrix(0,p,p)
# b[upper.tri(b)] <- runif(p*(p-1)/2, min = -b.magnitude, max = b.magnitude)
# b[abs(b) < lower.thresh] <- 0
# 
# zeros <- rbinom(p*(p-1)/2, 1, prob = prob_B)
# b[upper.tri(b)][which(zeros == 0)] <- 0
# s0 <- sum(abs(b)>1e-6)



# Given structure ---------------------------------------------------------
set.seed(394)
s0 <- sum(B_true != 0)
b <- matrix(0,p,p)
b[B_true != 0] <- runif(s0, min = -b.magnitude, max = b.magnitude)
dimnames(b) <- dimnames(B_true)





# generate Sigma (better way) ---------------------------------------

block.func <- function(x, eps){
  x[x == 1] <- runif(sum(x == 1), -eps, eps)
  diag(x) <- 1
  x[lower.tri(x)] <- 0
  x <- (x + t(x)) / 2
  inv.x <- round(solve(x), 10)
  inv.x <- round(cov2cor(inv.x), 10)
  x <- round(solve(inv.x), 10)
  return(x)
}

# switch(theta_type,
#        'random' = {
#          set.seed(89)
#          theta.val <- rbinom(n*n,1,prob = sparse.prob)
#          theta <- matrix(theta.val,n,n)
#          theta[theta == 1] <- runif(sum(theta == 1),-1,1)
#          diag(theta) <- 5
#          theta[lower.tri(theta)]<- 0
#          theta <- (theta+t(theta))/2
#          theta[abs(theta) < 0.05] = 0
#          if(!is.positive.definite(theta)) warning("Theta is not positive definite !!!")
#          sig <- round(solve(theta),10)
#          sig <- round(cov2cor(sig),10)
#          theta <- round(solve(sig),10)
#          theta[abs(theta) < 0.05] = 0
#        },
#        "block-diag" = {
#          set.seed(93)
#          m <- n%/%clique_size
#          m2 <- n%%clique_size
#          theta.val <- rbinom(clique_size^2*m + m2^2, 1, prob = sparse.prob*2)
#          diag.block <- list()
#          for(i in 1:m) diag.block[[i]] <- matrix(theta.val[(1 + clique_size^2*(i - 1)) : (clique_size^2*i)], clique_size, clique_size)
#          diag.block <- lapply(diag.block, block.func, eps = epsilon)
#          if (m2 != 0){
#            diag.block[[m + 1]] <- matrix(theta.val[(m*clique_size^2 + 1) : length(theta.val)], m2, m2)
#            diag.block[[m + 1]] <- block.func(diag.block[[m + 1]], eps = epsilon)
#          }
#          if(!all(unlist(lapply(diag.block, is.positive.definite))) == TRUE) warning("Theta is not positive definite!!")
#          theta <- bdiag(diag.block)
#          theta[abs(theta) < 0.05] = 0
#          sig <- round(solve(theta),10)
#          sig <- cov2cor(sig)
#          theta <- round(solve(sig), 10)
#        }
# )

set.seed(89)


# Toeplitz ----------------------------------------------------------------
# 
# first_block_size <- floor(n/2)
# second_block_size <- floor(n/4)
# sig <- diag(n)
# indices <- as.matrix(expand.grid(1:first_block_size, 1:first_block_size))
# sig.first <- diag(first_block_size)
# sig.second <- diag(second_block_size)
# for(i in 1:nrow(indices)) sig.first[indices[i,1], indices[i,2]] <- 0.3^(abs(diff(indices[i,]))/5)
# sig[1:first_block_size,1:first_block_size] <- sig.first
# indices2 <- as.matrix(expand.grid(1:second_block_size, 1:second_block_size))
# for(i in 1:nrow(indices2)) sig.second[indices2[i, 1], indices2[i,2]] <- 0.3^(abs(diff(indices2[i,]))/5)
# sig[(first_block_size + 1): (first_block_size + second_block_size), (first_block_size + 1): (first_block_size + second_block_size)] <- sig.second
# if(!is.positive.definite(sig)) warning("Theta is not positive definite !!!")
# theta <- round(solve(sig),10)
# theta[abs(theta) < 0.005] = 0


# Exp.Decay ---------------------------------------------------------------
# 
# first_block_size <- floor(n/3)
# second_block_size <- floor(n/2)
# theta <- diag(n)
# theta.first <- diag(first_block_size)
# theta.2nd <- diag(second_block_size)
# indices <- as.matrix(expand.grid(1:first_block_size, 1:first_block_size))
# for(i in 1:nrow(indices)) theta.first[indices[i,1], indices[i,2]] <- 0.3^(abs(diff(indices[i,]))/5)
# theta[1:first_block_size, 1:first_block_size] <- theta.first
# indices2 <- as.matrix(expand.grid(1:second_block_size, 1:second_block_size))
# for(i in 1:nrow(indices2)) theta.2nd[indices2[i,1], indices2[i,2]] <- 0.3^(abs(diff(indices2[i,]))/5)
# theta[first_block_size+1:second_block_size, first_block_size+1:second_block_size] <- theta.2nd
# sig <- round(cov2cor(solve(theta)),5)
# theta <- round(solve(sig), 5)

# Equi.cor ----------------------------------------------------------------
# 
# first_block_size <- floor(n/3)
# sec_block_size <- floor(n/4)
# sig <- diag(n)
# sig[1:first_block_size,1:first_block_size] <- 0.7
# sig[first_block_size+1:sec_block_size, first_block_size+1: sec_block_size] <- 0.7
# diag(sig) <- 1
# theta <- round(solve(sig),10)
# theta[abs(theta) < 0.005] = 0

Network811 <- read.table("real_Sigma/Cross_Parker-Consulting_info.txt", quote = "\"", comment.char = "")
freemans <- read.table("real_Sigma/Freemans_EIES-1_n48.txt", quote="\"", comment.char="")
USairport500 <- read.table("real_Sigma/USairport500.txt", quote="\"", comment.char="")

adjmat <- get.adjacency(graph_from_edgelist(as.matrix(freemans[,1:2]), directed = T))
n <- nrow(adjmat)
adjmat[adjmat == 2] <- 1
diag(adjmat) <- 1
sig <- matrix(0.9, n, n)
diag(sig) <- 1          
theta <- solve(sig)
theta[as.matrix(adjmat) != 1] <- 0
# sum(theta !=0)
sig <- cov2cor(solve(theta))
sig <- (t(sig) + sig)/2
theta <- round(solve(sig), 6)
sum(abs(theta) ==0)
zeropos <- which(adjmat ==0, arr.ind = T)
nonzero <- which(adjmat == 1, arr.ind = T)
set.seed(23)

# theta <- matrix(0, n, n)
# # val <- runif(nrow(zeropos), 0.1, 1) * sample(c(-1,1), size = nrow(zeropos), replace = T)
# for(i in 1:nrow(nonzero))  theta[nonzero[i, 1], nonzero[i, 2]] <- 1/(abs(nonzero[i, 1] - nonzero[i,2]) + 1)
# isSymmetric(theta)
is.positive.definite(as.matrix(theta))
# sig <- cov2cor(round(solve(theta), 5))
# sig <- (sig + t(sig)) /2
# theta <- solve(sig)
# theta[theta < 1e-4] <- 0



# 

if(!is.positive.definite(sig)) warning("Theta is not positive definite !!!")
# zeropos <- which(abs(as.matrix(theta)) < 0.001, arr.ind = T)




