main <- function(){
  rm(list = ls())
  # setwd("C:/Users/hangjian/Dropbox/research/code/")
  setwd("~/Dropbox/research/code")
  source("loadpackages.R")
  # helper functions --------------------------------------------------------
  source("helper_funcs.R")
  source("gen_params_funcs.R")
  source("summary_func.R")
  # models function ------------------------------------------------------------------
  source("lassoIdentTheta_func.R")
  source("main_iter_fast_func.R")
  source("scale_lasso_func.R")
  source("flipflop.R")
  # execution function ---------------------------------------------------------------
  source("run_three_model_func.R")
  source("run_many_sim_func.R")
  source("run_sim_main_func.R")
  source("BIC_flipflop_func.R")
  source('unordered_sims.R')
  source("main_iteration_fast2.R")
  source("ordered_sims.R")
  source("ordered_sims_new.R")
asdf

  # GES ---------------------------------------------------------------------
  args <- args_for_parameter()  
  # args$iid <- F
  dir.create(path = args$setting)
  setwd(args$setting)
  saveRDS(args, file = "args.rds")
  estimands <- generate_parameters(args = args, seed = 222, bname = "ecoli70",
                                   btype = "continuous", theta_name = NULL)
                                   # theta_name = "USairport500",
                                   # bname = "pathfinder", btype = "discrete")
  image(estimands$theta)
  image(as(estimands$b, class(estimands$theta)))
  estimands$sig <- diag(args$n)
  estimands$theta <- diag(args$n)
  # image(as(bestX$theta_est, class(estimands$theta)))
  saveRDS(estimands, file = "estimands.rds")
# run sims ----------------------------------------------------------------
  ordered_runsims_new(start = 2, 2, args, estimands, max.iter = 100, lam_div = 40)
  
  
  
  unordered_runsims(start_pos = 1, repp = 10, args, estimands, thr = 0.3,
                    max.iter = 20, lambda.len = 10, div = 120)
  
  
  
  settings <- c("072620", "072630", "072640", "072650", "072660")
  for(i in 1:length(settings)){
    add_sparsebn(sim_num = 10, setting = settings[i], today_date = as.Date("2019-7-28"), seed = 10)
  }
  add_sparsebn(start = 1, sim_num = 10, setting = "5035", today_date = as.Date("2019-5-4"), seed = 10)
  
}

# main()
