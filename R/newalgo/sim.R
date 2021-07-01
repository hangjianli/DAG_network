# load packages----------------------------------------------------------------
rm(list = ls())
setwd("~/Documents/research/dag_network/")
source("R/newalgo/loadpackages.R")
source("R/newalgo/helper_funcs.R")
source("R/newalgo/gen_params_funcs.R")
source("R/newalgo/functions.R")
source("R/newalgo/flipflop.R")


# generate specs ----------------------------------------------------------
args <- args_for_parameter()  
dir.create(path = paste0('output/', args$setting))
saveRDS(args, file = paste0('output/', args$setting, "/args.rds"))
estimands <- generate_parameters(args = args, seed = 5)
image(estimands$theta, axes=T)
image(as(estimands$b, class(estimands$theta)))
estimands$b[estimands$b!=0] %>% as.numeric() %>% hist(breaks=20)
estimands$theta[estimands$theta!=0] %>% as.numeric() %>% hist(breaks=20)
saveRDS(estimands, file = paste0('output/', args$setting, "/estimands.rds"))
estimands$s0
# X_ <- sim_X(
#   vers = 1,
#   n = args$n,
#   p = estimands$realp,
#   omg.sq = estimands$omg.sq,
#   sig = estimands$sig,
#   b = estimands$b
# )
# run simulation ----------------------------------------------------------
# ordered 
setwd("~/Documents/research/dag_network")
simID <- args$setting
sim_newalgo_ordered(
  args,
  estimands, 
  start_sim=3,
  end_sim=10, 
  lamLen=10,
  lambda2 = 20,
  lambda1_max_div = 10,
  lambda1_max_div2 = 400)

estimands$theta[1:20,1:20]
as(main_lam_6$thetahat[1:20,1:20], 'dgCMatrix')
image(estimands$theta, axes=T)

# simID = '121006'
# args <- readRDS(paste0("~/Documents/research/dag_network/output/", simID,"/args.rds"))
# estimands <- readRDS(paste0("~/Documents/research/dag_network/output/", simID, "/estimands.rds"))
process_output_ordered_tmp(
  simID=args$setting,
  estimands,
  args,
  start=5,
  num_sim=10,
  thr=0.1
)
get_average_shd_ordered(simID = args$setting, nsim = 10)


simID = args$setting
process_output_ordered(simID = args$setting, estimands = estimands, 
                       args = args, start = 3, num_sim = 10,  thr = 0.1)
get_all_shd_ordered(simID = args$setting, estimands, start = 3, num_sim = 10)
get_average_shd_ordered(simID = args$setting, nsim = 1)


get_roc_plots_noff(simID = args$setting, ymin = 50, xmax = 200, thrLen = 50, start = 1, nsim = 10)

# shd_average <- readRDS("~/Documents/research/dag_network/output/101/shd_average.rds")
# shd_average$shdXbaseline
# shd_average$shdXff
args <- readRDS(paste0("~/Documents/research/dag_network/output/", simID,"/args.rds"))
estimands <- readRDS(paste0("~/Documents/research/dag_network/output/", simID, "/estimands.rds"))
# get_roc_data_ordered(simID = simID)
get_all_shd_ordered(simID = simID, estimands, args$num_sim)
get_average_shd_ordered(simID = simID, nsim = as.numeric(args$num_sim))


shd_average <- readRDS("~/Documents/research/dag_network/output/5514/shd_average.rds")
xtable(t(shd_average), digits = c(0,1,1,1,3,3,1,1))



# unordered  --------------------------------------------------------------

rm(list = ls())
setwd("~/Documents/research/dag_network/")
source("R/loadpackages.R")
source("R/helper_funcs.R")
source("R/gen_params_funcs.R")
source("R/newalgo/functions.R")
source("R/newalgo/flipflop.R")

args <- args_for_parameter()  
dir.create(path = paste0('output/', args$setting))
saveRDS(args, file = paste0('output/', args$setting, "/args.rds"))
estimands <- generate_parameters(args = args, seed = 11)
image(estimands$theta, axes=T)
image(as(estimands$b, class(estimands$theta)))
estimands$b[estimands$b!=0] %>% as.numeric() %>% hist(breaks=20)
estimands$theta[estimands$theta!=0] %>% as.numeric() %>% hist(breaks=20)
saveRDS(estimands, file = paste0('output/', args$setting, "/estimands.rds"))
estimands$s0

estimands$theta[1:10,1:10]
image(as(main_lam_7$thetahat, class(estimands$theta)))
as(main_lam_7$thetahat[1:10,1:10], class(estimands$theta))

setwd("~/Documents/research/dag_network")
simID <- args$setting
sim_newalgo_unordered(
  args, 
  estimands, 
  start_sim=1, 
  end_sim=10,
  lamLen=10,
  lambda2 = 1,
  lambda1_max_div = 20,
  lambda1_max_div2 = 400)
process_output_unordered(simID = simID, thr = 0.1, start = 1, nsim = 10)
get_average_shd_unordered(simID = simID,start = 1, nsim = 10)


simIDarr = c(
  '00022', '00033', '00055', '00044',
  '321001', '321002', '321003', '321004', '321005'
)
simIDarr = c(
   # '5515','5505','5506', 
    '5508', '5514'
)

for(simID in simIDarr){
  args <- readRDS(paste0("~/Documents/research/dag_network/output/", simID,"/args.rds"))
  estimands <- readRDS(paste0("~/Documents/research/dag_network/output/", simID, "/estimands.rds"))
  process_output_unordered_tmp(
    simID = simID,
    estimands,
    args,
    start=1,
    nsim=args$num_sim, 
    thr = 0.1
  )
  get_average_shd_unordered(simID = simID, start = 1, nsim = 10)
  
}





simIDarr =c("5508", "5514")

shd_data <- get_shd_ji_diff(simIDarr)
# plot_boxplot_shd_ja(shd_data = shd_data)
shd_data

# library(xtable)
xtable(shd_average, digits = c(0,0,0))

mydata <- data.frame(
  shddiff = c(
    shd_data$shdgesdiff, 
    shd_data$shdpcdiff,
    shd_data$shdsbdiff
  ),
  jadiff = c(
    shd_data$JIgesdiff,
    shd_data$JIpcdiff,
    shd_data$JIsbdiff
  ),
  # id = rep(
  #   rep(
  #     rep(thetas, each=nsim), # 5 settings
  #     2), # (n < p) and (n ? p)
  #   3), # 3 methods
  label = rep(rep(c("n<p", "n>p"), each = 10), 3),
  method = rep(c("GES", "PC", "SBN"), each = nsim * 2)
)

setwd("~/Documents/research/dag_network/output/5508")

plot_shd_util(subset(mydata, (abs(shddiff) < 1000)), theta_type = 'andes_facebook',
              labels = c("n<p" = "(100, 192)", "n>p" = "(500, 192)"))
plot_ja_util(subset(mydata, (abs(jadiff) < 0.1)), theta_type = 'hepar2_celegan',
             labels = c("n<p" = "(100, 280)", "n>p" = "(500, 280)"))




# correct the theta errors for unordered sims -----------------------------

simIDs = c(
  '5501',
  '5502',
  '5503',
  '5515',
  '5505',
  '5506',
  '5508',
  '5514'
)


for(simID in simIDs){
  setwd(paste0('output/', simID))
  estimands <- readRDS('estimands.rds')
  support_size = sum(estimands$theta != 0)
  theta_err <- matrix(0, 10,2)
  dimnames(theta_err) = list(c(NULL, NULL), c('bcd', '1iter'))
  for(sim in 1:10){
    cat(paste0('[INFO] Processing sim ', sim, '\n'))
    bic_score <- readRDS(paste0(simID, '--', sim, '/BICscores_main.rds'))  
    bic_score_1iter <- readRDS(paste0(simID, '--', sim, '/BICscores_1iter_main.rds'))  
    best_bic <- which.min(bic_score)
    best_bic_1iter <- which.min(bic_score_1iter)
    best_res <- readRDS(paste0(simID, '--', sim, '/main_lam_', best_bic, '.rds'))
    best_res_1iter <- readRDS(paste0(simID, '--', sim, '/main_lam_', best_bic_1iter, '.rds'))  
    err1 = norm(estimands$theta - best_res$thetahat, type = '2')^2 / support_size
    err2 = norm(estimands$theta - best_res_1iter$thetahat, type = '2')^2 / support_size
    theta_err[sim,] = c(err1, err2)
  }
  saveRDS(theta_err, file = 'theta_err.rds')
  setwd("~/Documents/research/dag_network")
}




