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


generate_parameters <- function(args, seed = 10, bname = NULL, btype = NULL, 
                                theta_name = NULL){
  if(args$bchoice == "y") {
    cat("B was generated from data. \n")
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
    b.temp <- gen.B(p = args$p, seed = seed)
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
              # zerosb = zerosb, 
              theta = theta, sig = sig, 
              zeropos = zeropos, 
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
  omg <- runif(min = 0.1,max = 0.5,n = p) # omega's between 0,1
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
gen.B.from.btrue <- function(p, B_true, seed = 394, btype = NULL, lower.thresh = 0.1, b.mag = 1){
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

gen.theta.from.data <- function(name_theta, block_size = 20, n = 500, seed=5){
  mydata <- read.table(paste0("~/Dropbox/research/code/real_Sigma/", name_theta, ".txt"), quote="\"", comment.char="")
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
             sig <- cov2cor(solve(blocks[[i]]))
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
             sig <- cov2cor(solve(blocks[[i]]))
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
             blocks[[i]] <- matrix(1,block_size, block_size)
             diag(blocks[[i]])[2:block_size] <- sample(3:10, size = 1)
             sig <- cov2cor(blocks[[i]])
             blocks[[i]] <- round(solve(sig),5)
             
             if(!is.positive.definite(blocks[[i]])) stop("Theta is not positive definite !!!")
             zeropos[[i]] <- which(abs(as.matrix(blocks[[i]])) < 1e-3, arr.ind = T)
           }
           theta <- do.call(bdiag, blocks)
           theta[abs(theta) < 1e-3] = 0
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
             sig <- cov2cor(blocks[[i]])
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

