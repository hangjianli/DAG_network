res_pc <- rcausal::pc(df = XX, continuous = T,
                      depth = -1, significance = 0.01)

ptm <- proc.time()
res_pc <- rcausal::pc(df = XX, continuous = T, depth = -1, significance = 0.01)
# res_pc <- pcalg::pc(suffStat = list(C = cor(dataX), n = nrow(dataX)), indepTest = gaussCItest,
#                    alpha = 0.01, labels = var_name) 
run_time_pc <- proc.time() - ptm
# skl
adjmat_pc_CPDAG <- get_adjmat_from_pc(res_pc, args$p)
shd_pc <- compute_SHD_detail(adjmat_pc_CPDAG, adjmat_trueCPDAG)
