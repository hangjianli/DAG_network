for(tt in 1:nsim){
  setwd(dir = paste0(args$setting, "--", tt))
  today_date <- Sys.Date()-1
  
  BICscores_main <- readRDS(file = "BICscores_main.rds")
  bestk_bic_main <- which.min(BICscores_main)
  minrowcor_main <- readRDS(file = "minrowcor_main.rds")
  bestk_cor_main <- which.min(minrowcor_main)
  BICscores_bench <- readRDS(file = "BICscores_bench.rds")
  bestk_bic_bench <- which.min(BICscores_bench)
  minrowcor_bench <- readRDS(file = "minrowcor_bench.rds")
  bestk_cor_bench <- which.min(minrowcor_bench)
  BICscores_flipflop <- readRDS(file = "BICscores_flipflop.rds")
  bestk_bic_flip <- which.min(BICscores_flipflop)
  minrowcor_flip <- readRDS(file = "minrowcor_flip.rds")
  bestk_cor_flip <- which.min(minrowcor_flip)
  
  bestresult_flip_cor <- readRDS(file = paste0(format(today_date, "%Y-%m-%d"),
                                               "-nsim-", tt,
                                               "-vers-", bestk_cor_flip,
                                               "-lam-", bestk_cor_flip,
                                               "-flipflopMLE",".rds"))
  
  bestresult_flip_bic <- readRDS(file = paste0(format(today_date, "%Y-%m-%d"),
                                               "-nsim-", tt,
                                               "-vers-", bestk_bic_flip,
                                               "-lam-", bestk_bic_flip,
                                               "-flipflopMLE",".rds"))
  
  bestresult_main_bic <- readRDS(file = paste0(format(today_date, "%Y-%m-%d"),
                                               "-nsim-", tt,
                                               "-vers-", bestk_bic_main, "-lam-", bestk_bic_main,
                                               "-mainResult",".rds"))
  bestresult_bench_bic <- readRDS(file = paste0(format(today_date, "%Y-%m-%d"),
                                                "-nsim-", tt,
                                                "-vers-", bestk_bic_bench, "-lam-", bestk_bic_bench,
                                                "-benchMLEResult",".rds"))
  bestresult_main_cor <- readRDS(file = paste0(format(today_date, "%Y-%m-%d"),
                                               "-nsim-", tt,
                                               "-vers-", bestk_cor_main, "-lam-", 
                                               bestk_cor_main,
                                               "-mainMLEResult",".rds"))
  bestresult_bench_cor <- readRDS(file = paste0(format(today_date, "%Y-%m-%d"),
                                                "-nsim-", tt,
                                                "-vers-", bestk_cor_bench, "-lam-", 
                                                bestk_cor_bench,
                                                "-benchMLEResult",".rds"))
  
  thetaRE_mainbic <- norm(as.matrix(bestresult_main_bic$theta_est) - as.matrix(estimands$theta), "f") / norm(as.matrix(estimands$theta), "f")
  thetaRE_flipbic <- norm(bestresult_flip_bic$theta_est - estimands$theta, "f") / norm(estimands$theta, "f")
  
}