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