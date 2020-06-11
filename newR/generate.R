library(Matrix)
library(corpcor)
library(clusterGeneration)
# args_for_parameter <- function(){
#   # Input parameters from terminal ------------------------------------------
#   n <- readline(prompt = "Choose n >> ")
#   n <- as.integer(n)
#   p <- readline(prompt = "Choose p >> ")
#   p <- as.integer(p)
#   iid <- readline(prompt = "Equal variance for noise >> (T,F)  ")
#   iid <- as.logical(iid)
#   theta_type <- readline(prompt = "Choose theta structure: ('exp.decay', 'equi.cor', 'toeplitz', 'partially_connect', 'star', 'AR', 'random', 'diagonal') >> ")
#   while(!theta_type %in% c('exp.decay', 'diagonal', 'equi.cor', 'toeplitz','partially_connect', 'star', 'AR','random'))
#     theta_type <- readline(prompt = "Choose theta structure: ('exp.decay', 'equi.cor', 'toeplitz',
#                            'partially_connect', 'AR', 'random', 'diagonal') >> ")
#   block_size <- as.integer(readline(prompt = "Choose block size >> "))
#   bchoice <- readline(prompt = "Generate B from BN repository? (y or n) >> ")
#   while(!is.character(bchoice))
#     bchoice <- readline(prompt = "Generate B from BN repository? (y or n) >> ")
#   bchoice <- tolower(strsplit(bchoice, NULL)[[1]][1])
#   if(bchoice == "y"){
#     ncopy <- readline(prompt = "Choose number of copies of B >> ")
#     ncopy <- as.integer(ncopy)  
#   }else{
#     ncopy = 1
#   }
#   thetachoice <- readline(prompt = "Generate Theta from real networks? (y or n) >> ")
#   while(!is.character(thetachoice))
#     thetachoice <- readline(prompt = "Generate Theta from real networks? (y or n) >> ")
#   thetachoice <- tolower(strsplit(thetachoice, NULL)[[1]][1])
#   num_sim <- readline(prompt = "How many simulations do we run? >> ")
#   setting <- readline(prompt = "Assign a setting number (int). This is used to distingush different simulations >> ")
#   
#   return(list(n = n, 
#               p = p, 
#               ncopy = ncopy, 
#               iid = iid,
#               theta_type = theta_type,
#               block_size = block_size,
#               bchoice = bchoice,
#               thetachoice = thetachoice,
#               num_sim = num_sim,
#               setting = setting))
# }

genDAG <- function(p,)
#'     

genDAGdata <- function(args, 
                       seed = 10, 
                       bname = NULL, 
                       btype = NULL, 
                       theta_name = NULL
){
  if(args$bchoice == "y") {
    # cat("B was generated from data. \n")
    if(is.null(bname))
      bname = readline(prompt = "Pick a dataset: >> ")
    if(is.null(btype))
      btype = readline(prompt = "discrete or continuous? >> ")
    btrue <- try({constrB_from_BNrepo(name = bname, type = btype, ncopy = args$ncopy)})
    while(class(btrue) == 'try-error'){
      bname = readline(prompt = "Pick a dataset: >> ")
      btype = readline(prompt = "discrete or continuous? >> ")
      btrue <- try({constrB_from_BNrepo(name = bname, type = btype, ncopy = args$ncopy)})
    }
    pp = btrue$pp
    b.temp <- gen.B.from.btrue(p = btrue$realp, seed = seed, btype = btype,
                               B_true = btrue$B_true, lower.thresh = 0.1, b.mag = 1)
  }else{
    cat("B was generated from simulation.\n")
    b.temp <- gen.B(b.mag = 1, p = args$p, seed = seed)
    pp <- b.temp$pp
  }
  if(args$thetachoice == "y"){
    cat(paste0("Theta was generated from ", theta_name, " data set.\n"))
    theta.temp <- try({gen.theta.from.data(name_theta = theta_name, 
                                           block_size = args$block_size,
                                           n = args$n, seed = seed)})
  }else{
    cat("Theta was generated from simulation.\n")
    theta.temp <- gen.theta(struct = args$theta_type, 
                            n = args$n, 
                            block_size = args$block_size,
                            seed = seed*5)
  }
  
  # Must generate b first to determine p !!
  b <- b.temp$b
  s0 <- b.temp$s0
  # zerosb <- b.temp$zerosb
  realp <- b.temp$realp
  
  omg.temp <- gen.omg(realp, seed = seed, iid = args$iid)
  omg <- omg.temp$omg
  omg.sq <- omg.temp$omg.sq
  
  theta <- theta.temp$theta
  sig <- theta.temp$sig
  zeropos <- which(abs(as.matrix(theta)) < 1e-3, arr.ind = T)
  zeropos_list <- theta.temp$zeropos
  
  return(list(omg = omg, omg.sq = omg.sq, 
              b = b, s0 = s0,
              zerosb = zerosb,
              theta = theta, sig = sig, 
              zeropos = zeropos, 
              zeropos_list = zeropos_list,
              realp = realp, pp = pp,
              bname = bname, btype = btype,
              theta_name = theta_name))
}

# generate omega ----------------------------------------------------------
gen.omg <- function(p, iid=F, seed = 1){
  #' return the std of noise variable for all columns
  #' @param p, integer, number of nodes in DAG
  #' @param iid, boolean, default False
  #' @param seed, default 1
  #'
  
  if(iid) {
    return(rep(1,p))
  }
  set.seed(seed)
  omg <- sample(1:5, size = 10, replace = T) * 0.1
  return(omg)
}

# generate Beta-----------------------------------------------------

genB <- function(p,
                  nEdges,
                  ub = 1,
                  lb = 0.1,
                  permute = F,
                  seed = 1
                  ){
  #' generate random DAG B matrix with positive/negative entries
  #' @return matrix
  
  set.seed(seed)
  newDAG <- t(sparsebnUtils::random.dag(nnode = p, nedge = nEdges, permute = permute))
  newB <- (newDAG!=0) + 0
  newB[newB!=0] = runif(nEdges, lb, ub)*(2*rbinom(nEdges, 1, 0.5)-1)
  dimnames(newB) <- list(as.character(1:p), as.character(1:p))
  return(newB)
}


gen.B.from.btrue <- function(p,
                             B_true, 
                             seed = 394, 
                             btype = NULL, 
                             lower.thresh = 0.1, 
                             b.mag = 1){
  #' B_true is the adjacency matrix
  #' 
  if(btype == "continuous"){
    s0 <- sum(abs(B_true)>0)
    b <- B_true
    # dimnames(b) <- dimnames(B_true)
    return(list(b=b, 
                s0=s0,
                # zerosb=zerosb,
                realp=p))
  }
  set.seed(seed*3)
  # zerosb <- as.numeric(abs(as.numeric(B_true[upper.tri(B_true)])) > 1e-5)
  b <- matrix(0,p,p)
  s0 <- sum(abs(B_true)>0)
  mysigns <- 2*rbinom(s0,1,0.5) - 1
  b[B_true != 0] <- runif(s0, min = lower.thresh, max = b.mag) * mysigns
  # b[abs(b) < lower.thresh] <- 0
  # s0 <- sum(b!=0)
  dimnames(b) <- dimnames(B_true)
  return(list(b=b, 
              s0=s0,
              # zerosb=zerosb,
              realp=p))
}


constrB_from_BNrepo <- function(name = "andes", type = "discrete", ncopy = 1){
  # Given name and ncopy, generate adjacency matrix B_true
  
  load(paste0("~/Documents/research/dag_network/BNRepo/", name, ".rda"))
  # load("~/Dropbox/research/code/BNRepo/arth150.rda")
  # load(paste0("~/../Dropbox/research/code/BNRepo/", name, ".rda")) 
  ordering <- node.ordering(bn)
  pp <- length(ordering)
  
  if(type == "discrete"){
    B_true <- matrix(0, pp, pp)
    dimnames(B_true) <- list(ordering, ordering)
    for(i in 1:pp){
      if(length(bn[[ordering[i]]][[2]]) == 0) next
      B_true[bn[[ordering[i]]][[2]], ordering[i]] <- 1
    }
    if(ncopy > 1){
      Blist <- vector("list", ncopy)
      for(i in 1:ncopy) Blist[[i]] <- B_true
      B_true <- do.call(adiag, Blist)  
    }
    s0 <- sum(B_true == 1)
    p <- pp * ncopy  
    return(list(B_true = B_true, realp = p, pp = pp))
  }
  
  if(type == "continuous"){
    B_true <- matrix(0, pp, pp)
    dimnames(B_true) <- list(ordering, ordering)
    MAGIC_NIAB_B <- lapply(bn, "[[", 4)
    wrong_names <- names(MAGIC_NIAB_B)
    for(i in 1:length(MAGIC_NIAB_B)){
      if(length(names(MAGIC_NIAB_B[[i]])) == 1) next
      B_true[,wrong_names[i]][names(MAGIC_NIAB_B[[i]])[-1]] <- MAGIC_NIAB_B[[i]][-1]
    }
    
    if(ncopy > 1){
      Blist <- vector("list", ncopy)
      for(i in 1:ncopy) Blist[[i]] <- B_true
      B_true <- do.call(adiag, Blist)  
    }
    s0 <- sum(B_true[upper.tri(B_true)] != 0)
    p <- pp * ncopy
    return(list(B_true = B_true, realp = p, pp = pp))
  }
  stop("Bad network!!")
}




# generate Theta/Sigma ---------------------------------------
gen.theta.from.data <- function(name_theta, block_size = 20, n = 500, seed=5){
  mydata <- read.table(paste0("data/real_Sigma/", name_theta, ".txt"), quote="\"", comment.char="")
  adjmat <- get.adjacency(graph = graph_from_edgelist(el = as.matrix(mydata[,1:2]), directed = F))
  adjmat[adjmat == 2] <- 1
  diag(adjmat) <- 1
  trueN <- dim(adjmat)[1]
  if(trueN < n){
    cat("True N is ", trueN, "\n")
    stop("N is smaller than n from data. \n")
  }
  adjmat <- adjmat[1:n, 1:n]
  # image(as(adjmat, class(estimands$theta)))
  
  num_blocks = ceiling(n/block_size)
  blocks <- vector(mode = "list", length = num_blocks) # blocks of sigma
  zeropos <- vector(mode = "list", length = num_blocks)
  set.seed(seed)
  for(i in 1:num_blocks){
    samp_ind <- sample(n, block_size)
    # blocks[[i]] <- as.matrix(adjmat[(1 + (i-1)*block_size): min(i*block_size, n), (1 + (i-1)*block_size): min(i*block_size,n)])
    blocks[[i]] <- as.matrix(adjmat[samp_ind, samp_ind])
    # blocks[[i]][blocks[[i]] == 1] <- 0.7
    blocks[[i]][blocks[[i]] == 1] <- runif(n = sum(blocks[[i]]),-4,4)
    blocks[[i]] <- (blocks[[i]] + t(blocks[[i]])) / 2
    diag(blocks[[i]]) <- 1
    while(!is.positive.definite(blocks[[i]])) {
      diag(blocks[[i]]) <- diag(blocks[[i]]) + 1
      cat("Added 1 to diagonal. \n")
    }
    sig <- cov2cor(solve(blocks[[i]]))
    blocks[[i]] <- round(solve(sig),5)
    zeropos[[i]] <- which(abs(as.matrix(blocks[[i]])) < 1e-3, arr.ind = T)
  }
  theta <- do.call(bdiag, blocks)
  theta[abs(theta) < 1e-3] = 0
  sig <- round(solve(theta),5)
  
  return(list(theta=theta, sig=sig, theta_type=name_theta, zeropos=zeropos))
}

genToeplitz <- function(n,
                        nBlocks,
                        bSizes,
                        seed=1){
  #' generate AR(1) precision matrix
  #' @param n: integer, number of observations
  #' @param nBlocks: integer, number of blocks
  #' @param bSize: vector of size nBlocks, size of each block
  #' @return list: list of precision matrix and zero positions 
  blocks <- vector(mode = "list") 
  zeropos <- vector(mode = "list")
  
  set.seed(seed)
  for(i in 1:nBlocks){
    indices <- as.matrix(expand.grid(1:bSizes[i], 1:bSizes[i]))
    blocks[[i]] <- diag(bSizes[i])
    for(j in 1:nrow(indices)){
      blocks[[i]][indices[j,1], indices[j,2]] <- 0.7^abs(diff(indices[j,]))
    } 
    
    stopifnot(is.positive.definite(blocks[[i]])) 
    
    theta <- round(solve(blocks[[i]]),5)
    zeropos[[i]] <- which(abs(as.matrix(theta)) < 1e-3, arr.ind = T)
  }
  sig <- do.call(Matrix::bdiag, blocks)
  theta <- round(solve(sig),5)
  theta[abs(theta) < 1e-3] = 0
  return(list(theta=theta, zeropos=zeropos))
}


genExpDecay <- function(n,
                        nBlocks,
                        bSizes,
                        seed=1){
  #' generate exponential decay precision matrix
  
  blocks <- vector(mode = "list") 
  zeropos <- vector(mode = "list")
  set.seed(seed)
  for(i in 1:nBlocks){
    indices <- as.matrix(expand.grid(1:bSizes[i], 1:bSizes[i]))
    blocks[[i]] <- diag(bSizes[i])
    for(j in 1:nrow(indices)){
      blocks[[i]][indices[j,1], indices[j,2]] <- 0.7^abs(diff(indices[j,]))
    } 
    sig <- cov2cor(solve(blocks[[i]]))
    blocks[[i]] <- round(solve(sig),5)
    
    stopifnot(is.positive.definite(blocks[[i]])) 
    
    zeropos[[i]] <- which(abs(as.matrix(blocks[[i]])) < 1e-3, arr.ind = T)
  }
  theta <- do.call(bdiag, blocks)
  theta[abs(theta) < 1e-3] = 0
  sig <- round(solve(theta),5)
  
  return(list(theta=theta, zeropos=zeropos))
}


genRandomTheta <- function(n,
                           nBlocks,
                           bSizes,
                           prob=0.2,
                           seed=1){
  
  set.seed(seed)
  blocks <- vector(mode = "list") 
  zeropos <- vector(mode = "list")  
  
  for(i in 1:nBlocks){
    theta.val <- rbinom(bSizes[i]^2, 1, prob = prob)
    blocks[[i]] <- matrix(theta.val, bSizes[i], bSizes[i])
    blocks[[i]][blocks[[i]] == 1] <- runif(sum(blocks[[i]] == 1),-1,1)
    diag(blocks[[i]]) <- 5
    blocks[[i]][lower.tri(blocks[[i]])]<- 0
    blocks[[i]] <- (blocks[[i]]+t(blocks[[i]]))/2
    blocks[[i]][abs(blocks[[i]]) < 0.05] = 0
    
    sig <- cov2cor(solve(blocks[[i]]))
    blocks[[i]] <- round(solve(sig),5)
    
    stopifnot(is.positive.definite(blocks[[i]])) 
    
    zeropos[[i]] <- which(abs(blocks[[i]]) < 1e-3, arr.ind = T)
  }
  theta <- do.call(bdiag, blocks)
  theta[abs(theta) < 1e-3] = 0
  sig <- round(solve(theta),5)
  
  return(list(theta=theta, zeropos=zeropos))  
}

genAR <- function(n,
                  nBlocks,
                  bSizes,
                  bandWidthScale=4,
                  seed=1){
  
  blocks <- vector(mode = "list") 
  zeropos <- vector(mode = "list")
  set.seed(seed)
  for(i in 1:nBlocks){
    indices <- as.matrix(expand.grid(1:bSizes[i], 1:bSizes[i]))
    blocks[[i]] <- diag(bSizes[i])
    for(j in 1:nrow(indices)){blocks[[i]][indices[j,1], indices[j,2]] <- 0.7^(abs(diff(indices[j,])))} 
    blocks[[i]] <- band(blocks[[i]], k1 = -ceiling(bSizes[i]/bandWidthScale), k2 = ceiling(bSizes[i]/bandWidthScale))
    while(!is.positive.definite(as.matrix(blocks[[i]]))){
      cat(paste0("The ", i, "th block is not invertible! Adding 1. \n"))
      blocks[[i]] <- blocks[[i]] + 1*diag(dim(blocks[[i]])[1])
    }
    sig <- cov2cor(solve(blocks[[i]]))
    blocks[[i]] <- round(solve(sig),5)
    stopifnot(is.positive.definite(as.matrix(blocks[[i]]))) 
    zeropos[[i]] <- which(abs(as.matrix(blocks[[i]])) < 1e-3, arr.ind = T)
  }
  theta <- do.call(bdiag, blocks)
  theta[abs(theta) < 1e-3] = 0
  sig <- round(solve(theta),5)
  return(list(theta=theta, zeropos=zeropos))
}

genStar <- function(n,
                    nBlocks,
                    bSizes,
                    seed=1){
  #' generate star shaped precision matrix from a list of partial correlations
  #' 
  #' @return matrix: normalized to have correlation matrix as inverse
  blocks <- vector(mode = "list") # blocks of theta
  zeropos <- vector(mode = "list")
  set.seed(seed)
  
  for(i in 1:nBlocks){
    blocks[[i]] <- diag(bSizes[i])
    partial_cor <- runif(bSizes[i]-1, 0.1, 1) * sign(runif(1,-1,1))
    blocks[[i]][1,2:bSizes[i]] <- partial_cor
    blocks[[i]][2:bSizes[i], 1] <- partial_cor
    while(!is.positive.definite(as.matrix(blocks[[i]]))){
      blocks[[i]] <- blocks[[i]] + diag(bSizes[i])*.1
    }
    sig <- cov2cor(solve(blocks[[i]]))
    blocks[[i]] <- round(solve(sig),5)
    stopifnot(is.positive.definite(as.matrix(blocks[[i]]))) 
    zeropos[[i]] <- which(abs(as.matrix(blocks[[i]])) < 1e-3, arr.ind = T)
  }
  theta <- do.call(bdiag, blocks)
  theta[abs(theta) < 1e-3] = 0
  sig <- round(solve(theta),5)
  
  return(list(theta=theta, zeropos=zeropos))
}

genEquiCor <- function(n,
                       nBlocks,
                       bSizes,
                       corrs=c(0.6, 0.7),
                       seed=1){
  #' return an equi.cor covariance matrix with covariances 0.6,0.7
  #' @output a list 
  blocks <- vector(mode = "list") # blocks of theta
  zeropos <- vector(mode = "list")
  set.seed(seed)
  
  for(i in 1:nBlocks){
    blocks[[i]] <- matrix(sample(corrs,1)*sign(runif(1,-1,1)), bSizes[i], bSizes[i])
    diag(blocks[[i]]) <- 1
    while(!is.positive.definite(as.matrix(blocks[[i]]))){
      blocks[[i]] <- blocks[[i]] + diag(bSizes[i])*0.1
    }
    blocks[[i]] <- cov2cor(blocks[[i]])
    stopifnot(is.positive.definite(as.matrix(blocks[[i]]))) 
    theta <- round(solve(blocks[[i]]),5)
    zeropos[[i]] <- which(abs(as.matrix(theta)) < 1e-3, arr.ind = T)
  }
  sig <- do.call(bdiag, blocks)
  theta <- round(solve(sig),10)
  theta[abs(theta) < 1e-3] = 0
  
  return(list(theta=theta, zeropos=zeropos))
}


genPartiallyConnect <- function(n,
                                nBlocks,
                                bSizes,
                                threashold=0.05,
                                seed=1){
  #' generate partially connected precision matrix by threasholding entries from precision
  #' 
  blocks <- vector(mode = "list")
  zeropos <- vector(mode = "list")
  set.seed(seed)
  for(i in 1:nBlocks){
    blocks[[i]] <- genPositiveDefMat(covMethod = "eigen", dim = bSizes[i])$Sigma
    # blocks[[i]][abs(blocks[[i]]) < 0.1] <- 0
    sig <- cov2cor(blocks[[i]])
    blocks[[i]] <- round(solve(sig),5)
    blocks[[i]][abs(blocks[[i]]) < threashold] <- 0
    stopifnot(is.positive.definite(blocks[[i]]))
    zeropos[[i]] <- which(abs(as.matrix(blocks[[i]])) < 1e-3, arr.ind = T)
  }
  theta <- do.call(bdiag, blocks)
  theta[abs(theta) < 1e-3] = 0
  sig <- round(solve(theta),5)
  
  return(list(theta=theta, zeropos=zeropos))
}

genIdentity <- function(n,
                        nBlocks,
                        bSizes){
    blocks <- vector(mode = "list") 
    zeropos <- vector(mode = "list")
    for(i in 1:nBlocks){
      blocks[[i]] <- diag(bSizes[i])
      zeropos[[i]] <- which(abs(as.matrix(blocks[[i]])) < 1e-3, arr.ind = T)
    }
    sig <- do.call(bdiag, blocks)
    theta <- round(solve(sig),5)
    return(list(theta=theta, zeropos=zeropos))
}

genTheta <- function(n,
                     graphType = 'exp.decay',
                     prob = 1/5,
                     nBlocks = 5,
                     randomBS= T, 
                     bandScale = 4,
                     corrs = c(0.6,0.7),
                     threashold=0.05,
                     seed = 1
                     ){
  #' generate the precision matrix
  #' @param randomBS: boolean. Whether or not to use random block sizes.
  #' @param prob: float number, sparsity level in genRandomTheta
  #' @param bandScale:
  #' @param corrs: vector. correlations in genEquiCor
  #' @param threashold: float number, threashold for partially connect
  #' @return list of precision matrix and zero positions
  if(randomBS){
    set.seed(seed)
    prop <- runif(nBlocks) 
    block_size <- round(prop/sum(prop) * n, digits = 0)
    stopifnot(sum(block_size) == n)
  }else{
    block_size = rep(n / nBlocks, nBlocks)
  }
  num_blocks = ceiling(n/block_size)
  message('Block sizes: ', paste(block_size, collapse = " "))
  switch(graphType, 
         genIdentity = genIdentity(n, nBlocks, bSizes),
         genToeplitz = genToeplitz(n, nBlocks, bSizes, seed),
         genExpDecay = genExpDecay(n, nBlocks, bSizes, seed),
         genRandomTheta = genRandomTheta(n, nBlocks, bSizes, prob, seed),
         genAR = genAR(n, nBlocks, bSizes, bandScale, seed),
         genStar = genStar(n, nBlocks, bSizes, seed),
         genEquiCor = genEquiCor(n, nBlocks, bSizes, corrs, seed),
         genPartiallyConnect = genPartiallyConnect(n, nBlocks, bSizes, threashold, seed)
  )
}

# Simulate observational data from network DAG ----------------------------

sim_X <- function(vers, omg, sig, b){
  #' sample data from the given network DAG
  #' @param omg vector. Standard deviation of error
  #' @param sig matrix. row-covariance matrix 
  #' @param b matrix. DAG weights
  #'@return A list of X and E, both are matrices
  p = length(omg)
  n = dim(sig)[1]
  set.seed(vers)
  eps_mat <- matrix(0, n, p)
  eps_mat[,1] <- mvrnorm(1, mu = rep(0, n), Sigma = omg[1]^2*sig)
  eps_mat[,2] <- mvrnorm(1, mu = rep(0, n), Sigma = omg[2]^2*sig)
  X <- matrix(0, n, p)
  X[,1] <- eps_mat[,1]
  X[,2] <- X[,1]*b[1,2] + eps_mat[,2]
  for(i in 3:p) {
    eps_mat[, i] <- mvrnorm(1, mu = rep(0, n), Sigma = omg[i]^2*sig)
    X[,i] <- rowSums(sweep(X[,1:i-1], MARGIN = 2, b[1:i-1,i], "*")) + eps_mat[,i]
    if (i %% 50 == 0)
      cat("Getting ", i, "th column of X. \n" )
  }
  dimnames(X) <- list(NULL, as.character(1:p))
  return(list(X = X, eps_mat = eps_mat))
}


