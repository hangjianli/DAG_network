# load("BNRepo/magic-niab.rda") #44 continuous
# load("BNRepo/hailfinder.rda") #56
# load("BNRepo/win95pts.rda") #76
# load("BNRepo/hepar2.rda") #70
# load("BNRepo/barley.rda") #48
# load("BNRepo/insurance.rda") #27
# load("BNRepo/mildew.rda") #35
load("BNRepo/magic-irri.rda") #64 continuous
# load("BNRepo/arth150.rda") #107 continuous
# load("BNRepo/pathfinder.rda") #109
# load("BNRepo/ecoli70.rda") #46
# load("BNRepo/andes.rda") # 223
# load("BNRepo/diabetes.rda") #413



# Constructing B indicator matrix from bn.fit (continuous) ---------------------------------
ordering <- node.ordering(bn)
pp <- length(ordering)
MAGIC_NIAB_B <- lapply(bn, "[[", 4)
# # names(bn[[1]])
# # MAGIC_NIAB_B
parents <- lapply(bn, "[[", 2)
# root <- unlist(lapply(parents, length))[unlist(lapply(parents, length)) == 0]
# root
# 

# continuous network ------------------------------------------------------
B_true <- matrix(0, pp, pp)
dimnames(B_true) <- list(ordering, ordering)
wrong_names <- names(MAGIC_NIAB_B)
for(i in 1:length(MAGIC_NIAB_B)){
  if(length(names(MAGIC_NIAB_B[[i]])) == 1) next
  B_true[,wrong_names[i]][names(MAGIC_NIAB_B[[i]])[-1]] <- MAGIC_NIAB_B[[i]][-1]
}
s0 <- sum(B_true[upper.tri(B_true)] != 0)


# Constructing B indicator matrix from bn.fit (discrete) ---------------------------------
B_true <- matrix(0, pp, pp)
dimnames(B_true) <- list(ordering, ordering)
for(i in 1:pp){
  if(length(bn[[ordering[i]]][[2]]) == 0) next
  B_true[bn[[ordering[i]]][[2]], ordering[i]] <- 1
}

Blist <- vector("list", copies)
for(i in 1:copies) Blist[[i]] <- B_true
B_true <- do.call(adiag, Blist)
s0 <- sum(B_true == 1)
p <- pp * copies  
# # plot --------------------------------------------------------------------
# install.packages("Rgraphviz")
g <- Rgraphviz::layoutGraph(as.graphNEL(bn))
graph::nodeRenderInfo(g) <- list(fill="lightgreen", fontsize=20, 
                                 iwidth=list(a=7, d=3))
Rgraphviz::renderGraph(g)
