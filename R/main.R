rm(list = ls())
# setwd("~/Documents/research/dag_network/R")
source("R/loadpackages.R")
# helper functions --------------------------------------------------------
source("R/helper_funcs.R")
source("R/gen_params_funcs.R")
source("R/summary_func.R")
# models function ------------------------------------------------------------------
source("R/lassoIdentTheta_func.R")
source("R/main_iter_fast_func.R")
source("R/scale_lasso_func.R")
source("R/flipflop.R")
# execution function ---------------------------------------------------------------
source("R/run_three_model_func.R")
source("R/run_many_sim_func.R")
source("R/run_sim_main_func.R")
source("R/BIC_flipflop_func.R")
source('R/unordered_sims.R')
source("R/main_iteration_fast2.R")
source("R/ordered_sims.R")
source("R/ordered_sims_new.R")


p = 10
n = 5
omgTrue <- gen.omg(p = p, iid = T)
bTrue <- genB(p = p, nEdges = 1.5*p)
sigTrue <- genTheta


# GES ---------------------------------------------------------------------
args <- args_for_parameter()  
dir.create(path = paste0('output/', args$setting))
saveRDS(args, file = paste0('~/Documents/research/dag_network/output/', args$setting, "/args.rds"))

estimands <- generate_parameters(args = args, seed = 1, bname = "hailfinder",
                                 btype = "discrete", theta_name = "USairport500")
                                 # theta_name = "USairport500",
                                 # bname = "pathfinder", btype = "discrete")
image(estimands$theta)
image(as(estimands$b, class(estimands$theta)))
estimands$b[estimands$b!=0] %>% as.numeric() %>% hist(breaks=20)
# estimands$sig <- diag(args$n)
# estimands$theta <- diag(args$n)
# image(as(bestX$theta_est, class(estimands$theta)))
saveRDS(estimands, file = paste0('output/', args$setting, "/estimands.rds"))
# run sims ----------------------------------------------------------------
ordered_runsims_new(start = 2, 2, args, estimands, max.iter = 100, lam_div = 40)



unordered_runsims(start_pos = 1, repp = 10, args, estimands, thr = 0.3,
                  max.iter = 20, lambda.len = 10, div = 120)



settings <- c("072620", "072630", "072640", "072650", "072660")
for(i in 1:length(settings)){
  add_sparsebn(sim_num = 10, setting = settings[i], today_date = as.Date("2019-7-28"), seed = 10)
}
add_sparsebn(start = 1, sim_num = 10, setting = "5035", today_date = as.Date("2019-5-4"), seed = 10)
