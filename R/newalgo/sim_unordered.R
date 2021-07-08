# load packages----------------------------------------------------------------
rm(list = ls())
setwd("~/Documents/research/dag_network/")
source("R/newalgo/loadpackages.R")
source("R/newalgo/helper_funcs.R")
source("R/newalgo/gen_params_funcs.R")
source("R/newalgo/functions.R")
source("R/newalgo/flipflop.R")


# unordered  --------------------------------------------------------------

args <- args_for_parameter()  
dir.create(path = paste0('output/', args$setting))
saveRDS(args, file = paste0('output/', args$setting, "/args.rds"))
estimands <- generate_parameters(
  args = args, 
  # bname = 'hepar2',
  # btype = 'discrete',
  # theta_name = 'facebook',
  sbm_nblocks = 5,
  cluster_sizes = c(15, 40, 50, 80, 15),
  # cluster_sizes = c(20, 30, 5, 15,25, 5),
  sbm_probs = c(0.1,0.5),
  seed = 1)

  
image(estimands$theta, axes=T)
image(as(estimands$b, class(estimands$theta)))
estimands$b[estimands$b!=0] %>% as.numeric() %>% hist(breaks=20)
estimands$theta[estimands$theta!=0] %>% as.numeric() %>% hist(breaks=20)
saveRDS(estimands, file = paste0('output/', args$setting, "/estimands.rds"))
estimands$s0

setwd("~/Documents/research/dag_network")
simID <- args$setting 
sim_newalgo_unordered(
  args, 
  estimands, 
  start_sim=3, 
  end_sim=10,
  lamLen=10,
  lambda2 = 1,
  lambda1_max_div = 20,
  lambda1_max_div2 = 500)
process_output_unordered(simID = simID, thr = 0.1, start = 2, nsim = 10)
get_average_shd_unordered(simID = simID,start = 1, nsim = 10)


# for(simID in simIDarr){
#   args <- readRDS(paste0("~/Documents/research/dag_network/output/", simID,"/args.rds"))
#   estimands <- readRDS(paste0("~/Documents/research/dag_network/output/", simID, "/estimands.rds"))
#   process_output_unordered_tmp(
#     simID = simID,
#     estimands,
#     args,
#     start=1,
#     nsim=args$num_sim, 
#     thr = 0.1
#   )
#   get_average_shd_unordered(simID = simID, start = 1, nsim = 10)
# }





simIDarr =c("8004", "8005")

shd_data <- get_shd_ji_diff(simIDarr)
# plot_boxplot_shd_ja(shd_data = shd_data)
shd_data

# library(xtable)
# shd_average <- readRDS("~/Documents/research/dag_network/output/805/shd_average.rds")
# xtable(t(shd_average), digits = c(0,1,1,1,3,3,1,1,5))

nsim = 10
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
  method = rep(c("GES", "PC", "SBN"), each = 2*nsim)
)

setwd("~/Documents/research/dag_network/output/")

plot_shd_util(subset(mydata, (abs(shddiff) < 800)), theta_type = 'sbm',
              labels = c("n<p" = "(100, 300)", "n>p" = "(200, 100)"))
plot_ja_util(subset(mydata, (abs(jadiff) < 100)), theta_type = 'sbm',
             labels = c("n<p" = "(100, 300)", "n>p" = "(200, 100)"))


plot_ja_util2 <- function(df, theta_type="Equal corr", labels=c("n<p" = "(100, 200)", "n>p" = "(300, 100)")){
  scaleFUN <- function(x) sprintf("%.2f", x)
  setEPS()
  postscript(paste0(theta_type, "-ja.eps"))
  par(mar = c(1,1,1,1))
  plt <- ggplot(df, aes(y = jadiff)) + 
    geom_boxplot(aes(fill = method), width=1.5) + 
    # xlab("Sample size (n, p)") + 
    ylab("Increase in Jaccard index") + 
    # facet_grid(.~label) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 1.5) + 
    theme(
      # strip.text.x = element_text(size = 25),
      panel.background = element_blank(),
      panel.grid.minor = element_blank(),
      # panel.border = element_rect(color = "black", fill = NA),
      axis.line = element_line(colour = "black"),
      # axis.title=element_text(size=30,face="bold"),
      axis.title=element_blank(),
      axis.text = element_text(size=30),
      legend.text=element_text(size=30),
      axis.text.x=element_blank(),
      legend.title = element_text(size = 30),
      strip.background = element_rect(color = "black", fill = "white")
    ) + 
    scale_x_discrete(labels=labels) + 
    scale_y_continuous(labels = scaleFUN)
  print(plt)
  dev.off()
}

plot_shd_util2 <- function(
  df, 
  theta_type="Equal corr",
  labels=c("n<p" = "(100, 200)", "n>p" = "(300, 100)")
){
  scaleFUN <- function(x) sprintf("%.2f", x)
  setEPS()
  postscript(paste0(theta_type, "-shd.eps"))
  par(mar = c(1,1,1,1))
  plt <- ggplot(df, aes(y = shddiff)) + 
    geom_boxplot(aes(fill = method), width=1.5) + 
    # xlab("Sample size (n, p)") + 
    ylab("Decrease in SHD") + 
    # facet_grid(.~label) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 1.5) + 
    theme(
      # strip.text.x = element_text(size = 25),
      panel.background = element_blank(),
      panel.grid.minor = element_blank(),
      # panel.border = element_rect(color = "black", fill = NA),
      axis.line = element_line(colour = "black"),
      # axis.title=element_text(size=30,face="bold"),
      axis.title=element_blank(),
      axis.text = element_text(size=30),
      legend.text=element_text(size=30),
      axis.text.x=element_blank(),
      legend.title = element_text(size = 30),
      strip.background = element_rect(color = "black", fill = "white")
    )+
    scale_x_discrete(labels=labels) + 
    scale_y_continuous(labels = scaleFUN)
  print(plt)
  dev.off()
}

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




