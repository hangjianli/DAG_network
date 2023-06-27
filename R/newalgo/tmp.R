rm(list = ls())
setwd("~/Documents/research/dag_network/")
source("R/newalgo/loadpackages.R")
source("R/newalgo/helper_funcs.R")
source("R/newalgo/gen_params_funcs.R")
source("R/newalgo/functions.R")
source("R/newalgo/flipflop.R")
options(scipen = 7)


# n < p -------------------------------------------------------------------

# Rerun PC and GES with given order, assuming observations are independent.
# 101, 002, 003, 007

start_sim=1
end_sim=1

setwd("~/Documents/research/dag_network")
simID <- '101'
args <- readRDS(paste0("~/Documents/research/dag_network/output/", simID,"/args.rds"))
estimands <- readRDS(paste0("~/Documents/research/dag_network/output/", simID, "/estimands.rds"))
start_time = Sys.time()

for(sim in start_sim:end_sim){
  setwd(paste0("output/",args$setting, "/", args$setting, "--", sim))
  cat("------------------\n")
  cat(paste0("Start sim ", sim, " at ", round((Sys.time() - start_time) , 2), " seconds. \n"))
  # con <- file("test_jmlr.log")
  # sink(con, append=TRUE)
  # sink(con, append=TRUE, type="message")
  X_ = readRDS("X.rds")
  X <- X_$X
  res = bnlearn::pc.stable(
    x = as.data.frame(X),
    blacklist=ordering2blacklist(dimnames(X)[[2]]),
    alpha = 0.1,
  )
  cat(paste0("Encode structural constraints", " at ", round((Sys.time() - start_time) , 2), " minutes. \n"))
  fitted = bn.fit(x = res, data = as.data.frame(X), method = 'mle-g')
  saveRDS(fitted, file = paste0("pc_ordered.rds"))
  # estimate MLE for PC -----------------------------------------------------
  n = dim(X)[1]
  p = dim(X)[2]
  pc_res <- readRDS(file = 'pc_ordered.rds')
  pc_adj_b = matrix(0,p,p)
  for(i in 2:p){
    if(length(pc_res[[i]]$coefficients) > 1){
      pc_adj_b[as.integer(names(pc_res[[i]]$coefficients)[-1]), 
               i] = pc_res[[i]]$coefficients[-1]
    }
  }
  mle_result <- dag_mle_estimation(
    X = X,
    Bhat = pc_adj_b,
    Lhat = diag(n)
  )
  saveRDS(mle_result, file = paste0("pc_MLE_result.rds"))
  # BIC_result <- BIC_dag(
  #   X = X,
  #   block_idx = block_idx,
  #   bmle = mle_result$Bmle,
  #   omgmle = mle_result$omgmlesq,
  #   theta = diag(n)
  # )
  # saveRDS(BIC_1iter_result,  paste0('BIC_1iter_result_', k, '.rds'))
  # sink()
  # sink(type="message")
  cat(paste0("Experiment ", sim, " completed", " at ", round((Sys.time() - start_time) , 2), " minutes. \n"))
  setwd("~/Documents/research/dag_network")
}


simID = args$setting
process_output_ordered(
  simID = args$setting,
  estimands = estimands, 
  args = args, 
  start = 1, num_sim = 1,  thr = 0.1, 
  ff_flag = T,
  pc_flag = T
)

# 

get_all_shd_ordered(
  stats_file_name = 'SHDstats_JMLR.rds',
  simID = args$setting, 
  estimands,
  start = 1, 
  num_sim = 1,
  ff_flag = T, 
  pc_flag = T
)
get_average_shd_ordered(simID = args$setting, nsim = 1, ff_flag = T, pc_flag = T)



# library(bnlearn)
fitted[[100]]
plot(fitted)

# comparing with GES and PC -----------------------------------------------

# GES search in the space of equivalence classes, and there is at most one DAG compatible with a given order.
# We can compare 

# GES (121005)







# pc_sol(Xp, originalX = Xp, thetahat = diag(args$n), decor = F)

n <- dim(X_)[1]
p <- dim(X_)[2]
suffstat <- list(C = cor(X_), n = n)

res_pc <- pcalg::pc(
  suffStat = suffstat, 
  indepTest = gaussCItest, 
  alpha = 0.09,
  m.max = 5,
  labels = dimnames(X_)[[2]]
)
res_pc
plot(res_pc)

adjmat_pc_CPDAG <- as(res_pc, "amat")
dag_pc = pdag2dag(res_pc@graph)
dag_adj_pc = showAmat(dag_pc$graph)

dag_pc_a_diff_order = bnlearn::pdag2dag(as.bn(res_pc@graph), ordering = dimnames(X_)[[2]])
# dag_pc_a_diff_order2 = bnlearn::pdag2dag(as.bn(res_pc@graph), ordering = sample(dimnames(X_)[[2]]))



dag_pc_a_diff_order

dag_adj_pc = showAmat(dag_pc_a_diff_order$graph)

dag_pc_a_diff_order$arcs
pdag = getGraph(adjmat_fgesCPDAG_X)

# Adjacency Matrix G:
# G[i,j] = 1/2 if edge mark of edge i-j at j is head(kid)/tail(parent).
dag_adj = showAmat(dag$graph)
# compare -----------------------------------------------------------------


setwd("~/Documents/research/dag_network")
simID <- '121002'
args <- readRDS(paste0("~/Documents/research/dag_network/output/", simID,"/args.rds"))
estimands <- readRDS(paste0("~/Documents/research/dag_network/output/", simID, "/estimands.rds"))

bstar_adj <- 1*(abs(estimands$b) > 0)
res = list()
res2 = list()
for(sim in 1:10){
  X_ = readRDS(paste0("output/", simID, '/',simID, '--', sim, "/X.rds"))
  X = X_$X
  n <- dim(X)[1]
  p <- dim(X)[2]
  suffstat <- list(C = cor(X), n = n)
  res_pc <- pcalg::pc(
    suffStat = suffstat, 
    indepTest = gaussCItest, 
    alpha = 0.05,
    m.max = 5,
    labels = dimnames(X)[[2]]
  )
  
  adjmat_pc_CPDAG <- as(res_pc, "amat")
  dag_pc = pdag2dag(res_pc@graph,keepVstruct = T)
  dag_adj_pc = showAmat(dag_pc$graph)
  
  dag_adj_pc = bnlearn::pdag2dag(as.bn(res_pc@graph), ordering = dimnames(X)[[2]])
  # getGraph(dag_adj_pc)
  
  shd_pc <- compute_SHD_dag(adj1 = dag_adj_pc, adj_true = bstar_adj, estimands$s0) %>% unlist()
  res[[sim]] = shd_pc
  
  score <- new("GaussL0penObsScore", data=X, use.cpp=TRUE, lambda=12)
  ges.fit <- ges(score)
  # dag_ges = pdag2dag(ges.fit$essgraph)
  # ges.fit$repr$weight.mat()
  myadjmat <- matrix(0, p, p)
  dimnames(myadjmat) <- list(dimnames(X)[[2]], dimnames(X)[[2]])
  edge_set = ges.fit$essgraph$`.->.in.edges`
  for (i in 1:length(edge_set)){
    for (j in edge_set[[i]]){
      myadjmat[i, j] = 1
    }
  }
  
  pdag = getGraph(myadjmat)
  # dag_ges = pdag2dag(pdag)
  # dag_adj_ges = showAmat(dag_ges$graph)
  
  dag_adj_ges = bnlearn::pdag2dag(as.bn(pdag), ordering = dimnames(X)[[2]])
  shd_ges = compute_SHD_dag(adj1 = dag_adj_ges, adj_true = bstar_adj, estimands$s0) %>% unlist()
  res2[[sim]] = shd_ges
  
}


total <- data.frame(
  row.names = names(res[[1]]), 
  pc=rep(0, 6),
  ges=rep(0, 6)
)
for (sim in 1:10){
  total <- total + data.frame(pc=res[[sim]], ges=res2[[sim]])
}

total <- round(total / 10, 7)
total
# total$B0 = c(estimands$s0, rep(0,num_statistic-1))

