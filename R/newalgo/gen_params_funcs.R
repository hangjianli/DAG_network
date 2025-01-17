# Generates omega, B, Theta -----------------------------------------------
options(digits = 5)

args_for_parameter <- function(){
  # Input parameters from terminal ------------------------------------------
  n <- readline(prompt = "Choose n >> ")
  n <- as.integer(n)
  p <- readline(prompt = "Choose p >> ")
  p <- as.integer(p)
  iid <- readline(prompt = "Equal variance for noise >> (T,F)  ")
  iid <- as.logical(iid)
  theta_type <- readline(prompt = "Choose theta structure: ('exp.decay', 'equi.cor', 'toeplitz', 'partially_connect', 'star', 'AR', 'random', 'diagonal') >> ")
  while(!theta_type %in% c('exp.decay', 'diagonal', 'equi.cor', 'toeplitz','partially_connect', 'star', 'AR','random'))
    theta_type <- readline(prompt = "Choose theta structure: ('exp.decay', 'equi.cor', 'toeplitz',
                           'partially_connect', 'AR', 'random', 'diagonal') >> ")
  block_size <- as.integer(readline(prompt = "Choose block size >> "))
  bchoice <- readline(prompt = "Generate B from BN repository? (y or n) >> ")
  while(!is.character(bchoice))
    bchoice <- readline(prompt = "Generate B from BN repository? (y or n) >> ")
  bchoice <- tolower(strsplit(bchoice, NULL)[[1]][1])
  if(bchoice == "y"){
    ncopy <- readline(prompt = "Choose number of copies of B >> ")
    ncopy <- as.integer(ncopy)  
  }else{
    ncopy = 1
  }
  thetachoice <- readline(prompt = "Generate Theta from real networks? (y or n) >> ")
  while(!is.character(thetachoice))
    thetachoice <- readline(prompt = "Generate Theta from real networks? (y or n) >> ")
  thetachoice <- tolower(strsplit(thetachoice, NULL)[[1]][1])
  num_sim <- readline(prompt = "How many simulations do we run? >> ")
  setting <- readline(prompt = "Assign a setting number (int). This is used to distingush different simulations >> ")
  
  return(list(n = n, 
              p = p, 
              ncopy = ncopy, 
              iid = iid,
              theta_type = theta_type,
              block_size = block_size,
              bchoice = bchoice,
              thetachoice = thetachoice,
              num_sim = num_sim,
              setting = setting))
}

#' Load the settings for generating DAG parameters
#'  
load_settings <- function(){
  
}

generate_parameters <- function(
  args, 
  seed = 10, 
  bname = NULL, 
  btype = NULL, 
  sbm_probs = c(0.1,0.5),
  sbm_nblocks = 5,
  cluster_sizes,
  theta_name = NULL
){
  if(args$bchoice == "y") {
    cat("B will be generated from data... \n")
    if(is.null(bname))
      bname = readline(prompt = "Pick a dataset from bnlearn: >> ")
    if(is.null(btype))
      btype = readline(prompt = "Is B discrete or continuous? >> ")
    btrue <- try({constrB_from_BNrepo(name = bname, type = btype, ncopy = args$ncopy)})
    while(class(btrue) == 'try-error'){
      bname = readline(prompt = "Pick a dataset from bnlearn: >> ")
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
    cat("Theta will be generated from data... \n")
    if(is.null(theta_name))
      theta_name = readline(prompt = "Pick a dataset for theta: >> ")
    cat(paste0("Theta was generated from ", theta_name, " data set.\n"))
    theta.temp <- try({gen.theta.from.data(
      name_theta = theta_name,
      theta_sparsity =theta_sparsity,
      block_size = args$block_size,
      n = args$n, seed = seed)})

  }else{
    cat("Theta was generated from simulation.\n")
    theta.temp <- gen.theta(struct = args$theta_type,
                            n = args$n,
                            block_size = args$block_size,
                            seed = seed*5)
    # theta.temp <- try({generate_sbm_theta(
    #   n = args$n,
    #   nblocks = sbm_nblocks,
    #   cluster_sizes = cluster_sizes,
    #   probs = sbm_probs,
    #   seed = seed)})
  }
  
  # Must generate b first to determine p !!
  b <- b.temp$b
  s0 <- b.temp$s0
  zerosb <- b.temp$zerosb
  realp <- b.temp$realp
  
  omg.temp <- gen.omg(realp, seed = seed, iid = args$iid)
  omg <- omg.temp$omg
  omg.sq <- omg.temp$omg.sq
  
  theta <- theta.temp$theta
  sig <- theta.temp$sig
  # zeropos <- which(abs(as.matrix(theta)) < 1e-3, arr.ind = T)
  zeropos_list <- theta.temp$zeropos
  
  return(list(omg = omg, omg.sq = omg.sq, 
              b = b, s0 = s0,
              zerosb = zerosb,
              theta = theta, sig = sig, 
              # zeropos = zeropos, 
              zeropos_list = zeropos_list,
              realp = realp, pp = pp,
              bname = bname, btype = btype,
              theta_name = theta_name))
}

# generate omega ----------------------------------------------------------
gen.omg <- function(p, seed = 23, iid=FALSE){
  if(iid){
    omg <- omg.sq <- rep(1,p)
    return(list(omg = omg, omg.sq = omg.sq))
  }
  set.seed(seed)
  omg <- runif(min = 0.3,max = 2,n = p) # omega's between 0,1
  omg[1] <- 1  
  omg.sq <- omg^2
  return(list(omg = omg, omg.sq = omg.sq))
}

# generate Beta-----------------------------------------------------

gen.B <- function(p, b.mag = 1, s0= 2*p, seed = 482, lower.thresh = 0.1){
  if(p <= 5 ) stop("p is too small!")
  set.seed(seed*2)
  bb <- randomDag(seed = seed, numNodes = p,numEdges = s0)
  b <-  get_adjmat_from_fges(bb$edges, length(bb$nodes), bb$nodes)
  b[b!=0] = runif(length(bb$edges), lower.thresh, b.mag)*(2*rbinom(length(bb$edges),1,0.5)-1)
  realp <- pp <- p
  dimnames(b) <- list(as.character(1:realp), as.character(1:realp))
  return(list(b=b, s0=s0, realp = realp, pp = pp))
}


## Given structure 
gen.B.from.btrue <- function(p, B_true, seed = 394, btype = NULL, lower.thresh = 0.5, b.mag = 0.8){
  # B_true is the adjacency matrix
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
  #' Given name and ncopy, generate adjacency matrix B_true
  #' 

  load(paste0("data/BNRepo/", name, ".rda"))
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


generate_sbm_theta <- function(n, nblocks, cluster_sizes, probs=c(0.1, 0.3), seed=1){
  cat("Generating theta from stochastic block model.. \n")
  P =  matrix(probs[1],n,n)
  # for(i in 1:nblocks){
  #   P[(1 + (i-1) * (ceiling(n/nblocks))): min(i * (ceiling(n/nblocks)), n), 
  #     (1 + (i-1) * (ceiling(n/nblocks))): min(i * (ceiling(n/nblocks)), n)] = probs[2]
  # }
  for(i in 1:length(cluster_sizes)){
    if(i==1){
      pre_sz = 0
    }else{
      pre_sz = pre_sz + cluster_sizes[i-1]
    }
    sz = cluster_sizes[i]
    P[(1+pre_sz):(pre_sz + sz),(1+pre_sz):(pre_sz + sz)] = probs[2]
  }
  
  diag(P) = 0
  set.seed(seed)
  A = matrix(rbinom(matrix(1,n,n),matrix(1,n,n),P),n,n)
  A = ceiling((A + t(A)) / 2)
  theta = diag(n)
  theta[A == 1] <- runif(n = sum(A),-5,5)
  
  # indices <- as.matrix(expand.grid(1:n, 1:n))
  # for(j in 1:nrow(indices)){
  #   theta[indices[j,1], indices[j,2]] <- 0.7^(abs(diff(indices[j,]))/5)
  # }
  # theta[A != 0] = 0
  theta = (theta + t(theta)) / 2
  theta = theta - (min(eigen(theta)$value)-0.1) * diag(dim(theta)[1])
  theta = as(theta, 'dgCMatrix')
  image(theta)
  
  # image(as(theta, 'dgCMatrix'))
  theta = round(theta, 5)  
  zeropos <- vector(mode = "list")
  zeropos[[1]] <- which(abs(as.matrix(theta)) < 1e-5, arr.ind = T)
  sig <- round(solve(theta),5)
  if(!is.positive.definite(as.matrix(sig))) stop("Sigma is not positive definite !!!")
  
  return(list(theta=theta, sig=sig, zeropos=zeropos))
}



gen.theta.from.data <- function(name_theta, theta_sparsity=0.7, block_size = 20, n = 500, seed=5){
  mydata <- read.table(paste0("data/real_Sigma/", name_theta, ".txt"), quote="\"", comment.char="")
  adjmat <- get.adjacency(graph = graph_from_edgelist(el = as.matrix(mydata[,1:2]), directed = F))
  adjmat[adjmat == 2] <- 1
  diag(adjmat) <- 1
  N <- dim(adjmat)[1]
  # set.seed(seed)
  # start = sample(N-n, 1)
  # adjmat <- adjmat[start:(start+n-1), start:(start+n-1)]
  # image(as(adjmat, class(estimands$theta)))
  
  num_blocks = ceiling(n/block_size)
  blocks <- vector(mode = "list", length = num_blocks) # blocks of sigma
  zeropos <- vector(mode = "list", length = num_blocks)
  set.seed(seed)
  
  # cursize = n
  
  for(i in 1:num_blocks){
    # block_size = min(cursize, block_size)
    # cursize = cursize - block_size
    samp_ind <- sample(1:N, block_size, replace = F)
    # theta <- as.matrix(adjmat[(1 + (i-1)*block_size): min(i*block_size, n), (1 + (i-1)*block_size): min(i*block_size,n)])
    theta <- as.matrix(adjmat[samp_ind, samp_ind])
    # blocks[[i]][blocks[[i]] == 1] <- 0.7
    theta_tmp <- diag(block_size)
    indices <- as.matrix(expand.grid(1:block_size, 1:block_size))
    for(j in 1:nrow(indices)) theta_tmp[indices[j,1], indices[j,2]] <- theta_sparsity^(abs(diff(indices[j,]))/3)
    theta_tmp[theta == 0] = 0
    theta = (theta_tmp + t(theta_tmp))/2
    # theta[theta == 1] <- runif(n = sum(theta),-5,5)
    # theta <- (theta+ t(theta)) / 2
    theta = theta - (min(eigen(theta)$value)-0.1) * diag(dim(theta)[1])
    # image(as(theta, class(estimands$theta)))
    # sig <- cov2cor(solve(theta))
    # blocks[[i]] <- round(solve(sig),5)
    blocks[[i]] <- round(solve(theta), 3)
    zeropos[[i]] <- which(abs(as.matrix(blocks[[i]])) < 1e-3, arr.ind = T)
  }
  theta <- do.call(bdiag, blocks)
  theta <- theta[1:n, 1:n]
  # theta[abs(theta) < 1e-4] = 0
  sig <- round(solve(theta),3)
  if(!is.positive.definite(as.matrix(sig))) 
    stop("Sigma is not positive definite !!!")
  return(list(theta=theta, sig=sig, theta_type=name_theta, zeropos=zeropos))
}


gene_theta_from_data2 <- function(
  
){
  mydata <- read.table(paste0("data/real_Sigma/", name_theta, ".txt"), quote="\"", comment.char="")
  adjmat <- get.adjacency(graph = graph_from_edgelist(el = as.matrix(mydata[,1:2]), directed = F))
  adjmat[adjmat == 2] <- 1
  diag(adjmat) <- 1
  N <- dim(adjmat)[1]
}


gen.theta <- function(struct = 'exp.decay',
                      n, 
                      sparse.prob = 1/5, 
                      # num_blocks = 3,
                      block_size = 10, #
                      fix.zero = T, 
                      seed = 89,
                      epsilon = 1){
  
  num_blocks = ceiling(n/block_size)
  
  switch(struct, 
         "diagonal" ={
           blocks <- vector(mode = "list") # blocks of sigma
           zeropos <- vector(mode = "list")
           for(i in 1:num_blocks){
             blocks[[i]] <- diag(block_size)
             zeropos[[i]] <- which(abs(as.matrix(blocks[[i]])) < 1e-3, arr.ind = T)
           }
           sig <- do.call(bdiag, blocks)
           theta <- round(solve(sig),5)
         }, 
         "toeplitz" = {
           blocks <- vector(mode = "list") # blocks of sigma
           zeropos <- vector(mode = "list")
           set.seed(seed)
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
           theta <- round(solve(sig),5)
           theta[abs(theta) < 1e-3] = 0
         },
         "exp.decay" = {
           blocks <- vector(mode = "list") # blocks of theta
           zeropos <- vector(mode = "list")
           set.seed(seed)
           for(i in 1:num_blocks){
             if(i == num_blocks)
               block_size = n - (num_blocks-1)*block_size
             indices <- as.matrix(expand.grid(1:block_size, 1:block_size))
             blocks[[i]] <- diag(block_size)
             for(j in 1:nrow(indices)) blocks[[i]][indices[j,1], indices[j,2]] <- 0.7^(abs(diff(indices[j,]))/5)
             # sig <- cov2cor(solve(blocks[[i]]))
             sig <- solve(blocks[[i]])
             blocks[[i]] <- round(solve(sig),5)
             if(!is.positive.definite(blocks[[i]])) stop("Theta is not positive definite !!!")
             zeropos[[i]] <- which(abs(as.matrix(blocks[[i]])) < 1e-3, arr.ind = T)
           }
           theta <- do.call(bdiag, blocks)
           theta[abs(theta) < 1e-3] = 0
           sig <- round(solve(theta),5)
         }, 
         "AR" = {
           blocks <- vector(mode = "list") # blocks of theta
           zeropos <- vector(mode = "list")
           set.seed(seed)
           for(i in 1:num_blocks){
             if(i == num_blocks)
               block_size = n - (num_blocks-1)*block_size
             indices <- as.matrix(expand.grid(1:block_size, 1:block_size))
             blocks[[i]] <- diag(block_size)
             for(j in 1:nrow(indices)) blocks[[i]][indices[j,1], indices[j,2]] <- 0.7^(abs(diff(indices[j,])))
             blocks[[i]] <- band(blocks[[i]], k1 = -ceiling(block_size/4), k2 = ceiling(block_size/4))
             while(!is.positive.definite(as.matrix(blocks[[i]]))){
               cat(paste0("The ", i, "th block is not invertible! Adding 1. \n"))
               blocks[[i]] <- blocks[[i]] + 1*diag(dim(blocks[[i]])[1])
             }
             # sig <- cov2cor(solve(blocks[[i]]))
             sig <- solve(blocks[[i]])
             blocks[[i]] <- round(solve(sig),5)
             if(!is.positive.definite(as.matrix(blocks[[i]]))) stop("Theta is not positive definite !!!")
             zeropos[[i]] <- which(abs(as.matrix(blocks[[i]])) < 1e-3, arr.ind = T)
           }
           theta <- do.call(bdiag, blocks)
           theta[abs(theta) < 1e-3] = 0
           sig <- round(solve(theta),5)
         },
         "star" = {
           blocks <- vector(mode = "list") # blocks of theta
           zeropos <- vector(mode = "list")
           set.seed(seed)
           for(i in 1:num_blocks){
             if(i == num_blocks)
               block_size = n - (num_blocks-1)*block_size
             blocks[[i]] <- matrix(0,block_size, block_size)
             diag(blocks[[i]]) = 1
             blocks[[i]][1, 2:block_size] = (1 / (block_size-1)) - 0.01
             blocks[[i]][2:block_size, 1] = (1 / (block_size-1)) - 0.01
             theta <- blocks[[i]]
             if(!is.positive.definite(blocks[[i]])) stop("Theta is not positive definite !!!")
             zeropos[[i]] <- which(abs(as.matrix(blocks[[i]])) < 1e-3, arr.ind = T)
           }
           theta <- do.call(bdiag, blocks)
           sig <- round(solve(theta),5)
         },
         "equi.cor" = {
           blocks <- vector(mode = "list") # blocks of theta
           zeropos <- vector(mode = "list")
           set.seed(seed)
           for(i in 1:num_blocks){
             if(i == num_blocks)
               block_size = n - (num_blocks-1)*block_size
             blocks[[i]] <- matrix(0.7,block_size, block_size)
             diag(blocks[[i]]) <- 1
             if(!is.positive.definite(blocks[[i]])) stop("Theta is not positive definite !!!")
             zeropos[[i]] <- which(abs(as.matrix(blocks[[i]])) < 1e-3, arr.ind = T)
             print(zeropos[[i]])
           }
           sig <- do.call(bdiag, blocks)
           theta <- round(solve(sig),10)
           theta[abs(theta) < 1e-3] = 0
         },
         "partially_connect" = {
           blocks <- vector(mode = "list") # blocks of theta
           zeropos <- vector(mode = "list")
           set.seed(seed)
           for(i in 1:num_blocks){
             if(i == num_blocks)
               block_size = n - (num_blocks-1)*block_size
             blocks[[i]] <- genPositiveDefMat(covMethod = "eigen", dim = block_size)$Sigma
             # blocks[[i]][abs(blocks[[i]]) < 0.1] <- 0
             # sig <- cov2cor(blocks[[i]])
             sig <- blocks[[i]]
             blocks[[i]] <- round(solve(sig),5)
             blocks[[i]][abs(blocks[[i]]) < 0.1] <- 0
             if(!is.positive.definite(blocks[[i]]))
               stop(paste0("Block ",i," is not PD!"))
             zeropos[[i]] <- which(abs(as.matrix(blocks[[i]])) < 1e-3, arr.ind = T)
           }
           theta <- do.call(bdiag, blocks)
           theta[abs(theta) < 1e-3] = 0
           sig <- round(solve(theta),5)
         },
         "block_diag_dense" = {
           blocks <- vector(mode = "list") # blocks of theta
           zeropos <- vector(mode = "list")
           set.seed(seed)
           for(i in 1:num_blocks){
             if(i == num_blocks)
               block_size = n - (num_blocks-1)*block_size
             A <- matrix(runif(block_size^2)*2-1, ncol=block_size) 
             blocks[[i]] <- t(A) %*% A
             # sig <- cov2cor(blocks[[i]])
             sig <- blocks[[i]]
             blocks[[i]] <- round(solve(sig),5)
             blocks[[i]][abs(blocks[[i]]) < 0.1] <- 0
             if(!is.positive.definite(blocks[[i]]))
               stop(paste0("Block ",i," is not PD!"))
             zeropos[[i]] <- which(abs(as.matrix(blocks[[i]])) < 1e-3, arr.ind = T)
           }
           theta <- do.call(bdiag, blocks)
           theta[abs(theta) < 1e-3] = 0
           sig <- round(solve(theta),5)
         },
         "random" = {
           set.seed(seed)
           # theta.val <- rbinom(n*n,1,prob = sparse.prob)
           # theta <- matrix(theta.val,n,n)
           
           theta.val1 <- rbinom(first_block_size^2,1,prob = sparse.prob)
           theta.val2 <- rbinom(second_block_size^2,1, prob = sparse.prob)
           firstblock <- matrix(theta.val1, first_block_size, first_block_size)
           secondblock <- matrix(theta.val2, second_block_size, second_block_size)
           theta <- matrix(0,n,n)
           theta[1:first_block_size, 1:first_block_size] <- firstblock
           theta[first_block_size+1:second_block_size, first_block_size+1:second_block_size] <- secondblock
           theta[theta == 1] <- runif(sum(theta == 1),-1,1)
           diag(theta) <- 5
           theta[lower.tri(theta)]<- 0
           theta <- (theta+t(theta))/2
           theta[abs(theta) < 0.05] = 0
           if(!is.positive.definite(theta)) warning("Theta is not positive definite !!!")
           sig <- round(solve(theta),10)
           sig <- round(cov2cor(sig),10)
           theta <- round(solve(sig),10)
           theta[abs(theta) < 0.05] = 0
           }
         
  )
  return(list(theta=theta, sig=sig, theta_type=struct, zeropos=zeropos))
}




# 
# # Network811 <- read.table("real_Sigma/Cross_Parker-Consulting_info.txt", quote = "\"", comment.char = "")
# freemans <- read.table("real_Sigma/Freemans_EIES-1_n48.txt", quote="\"", comment.char="")
# # USairport500 <- read.table("real_Sigma/USairport500.txt", quote="\"", comment.char="")


gen.theta.from.adj <- function(theta_name){
  raw_theta = read.table(paste0("real_Sigma/" , theta_name, ".txt"), quote = "\"", comment.char = "")
  adjmat <- get.adjacency(graph_from_edgelist(as.matrix(raw_theta[,1:2]), directed = F))
  adjmat[adjmat!=0]=1
  if(!isSymmetric(adjmat)){
    stop("Theta is not symmetric!")
  }
  
  n <- nrow(adjmat)
  diag(adjmat) <- 1
  sig <- matrix(0.9, n, n)
  # diag(sig) <- 1          
  # theta <- solve(sig)
  # theta[as.matrix(adjmat) != 1] <- 0
  # # sum(theta !=0)
  # sig <- cov2cor(solve(theta))
  # sig <- (t(sig) + sig)/2
  # theta <- round(solve(sig), 6)
  # sum(abs(theta) ==0)
  # zeropos <- which(adjmat ==0, arr.ind = T)
  # nonzero <- which(adjmat == 1, arr.ind = T)
  return(theta)
}

# adjmat <- get.adjacency(graph_from_edgelist(as.matrix(USairport500[,1:2]), directed = F))
# n <- nrow(adjmat)
# adjmat[adjmat == 2] <- 1
# diag(adjmat) <- 1
# sig <- matrix(0.9, n, n)
# diag(sig) <- 1          
# theta <- solve(sig)
# theta[as.matrix(adjmat) != 1] <- 0
# # sum(theta !=0)
# sig <- cov2cor(solve(theta))
# sig <- (t(sig) + sig)/2
# theta <- round(solve(sig), 6)
# sum(abs(theta) ==0)
# zeropos <- which(adjmat ==0, arr.ind = T)
# nonzero <- which(adjmat == 1, arr.ind = T)
# set.seed(23)

# theta <- matrix(0, n, n)
# # val <- runif(nrow(zeropos), 0.1, 1) * sample(c(-1,1), size = nrow(zeropos), replace = T)
# for(i in 1:nrow(nonzero))  theta[nonzero[i, 1], nonzero[i, 2]] <- 1/(abs(nonzero[i, 1] - nonzero[i,2]) + 1)
# isSymmetric(theta)
# sig <- cov2cor(round(solve(theta), 5))
# sig <- (sig + t(sig)) /2
# theta <- solve(sig)
# theta[theta < 1e-4] <- 0



# 

# if(!is.positive.definite(sig)) warning("Theta is not positive definite !!!")
# zeropos <- which(abs(as.matrix(theta)) < 0.001, arr.ind = T)



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

