library(ggplot2)


dat <- data.frame(Beta=c("andes", "hailfinder", "barley"),
                  Sigma=c("facebook", "celegans_n306"), 
                  method=c())

# ordered -----------------------------------------------------------------
datmain <- data.frame(Beta=rep(c("andes", "hailfinder", "barley"),2),
                     Sigma=rep(c("facebook", "celegans_n306"),each=3),
                     SHD = c(215.5, 204.30, 276.60, 260.80,213.50, ))

datbench <- data.frame(Beta=rep(c("andes", "hailfinder", "barley"),2),
                      Sigma=rep(c("facebook", "celegans_n306"),each=3),
                      SHD = c())


num_sim <- args$num_sim
threshes <- seq(0, 0.8,length.out = 20)
allShdS <- vector(mode = "list", length = length(threshes))

# today_date <- as.Date("2019-05-02")
today_date <- Sys.Date()


for(sim in 1:1){
  setwd(paste0("~/Dropbox/research/code/", args$setting,"--", sim))
  estimands <- readRDS(paste0("~/Dropbox/research/code/", args$setting,"--", sim, "/estimands.rds"))
  X_ <- readRDS(file = paste0(format(today_date, "%Y-%m-%d"), "-vers-",
                            args$setting, "-X-",sim,".rds"))
  X <- X_$X
  n <- dim(X)[1]
  p <- dim(X)[2]
  # dimnames(estimands$b) <- list(as.character(1:p), as.character(1:p))
  # adjmat_true_dag <- 1*(estimands$b != 0)
  # saveRDS(adjmat_true_dag, "adjmat_true_dag.rds")
  adjmat_trueCPDAG <- readRDS(paste0("~/Dropbox/research/code/", args$setting,"--", sim, "/adjmat_trueCPDAG.rds"))
  
  # BICscores_bench <- readRDS(paste0("~/Dropbox/research/code/",args$setting, "--",sim,"/BICscores_bench.rds"))
  BICscores_main <- readRDS(paste0("~/Dropbox/research/code/",args$setting, "--", sim,"/BICscores_main.rds"))
  minrowcor_main <- readRDS(paste0("~/Dropbox/research/code/",args$setting, "--", sim,"/minrowcor_main.rds"))
  # minrowcor_bench <- readRDS(paste0("~/Dropbox/research/code/",args$setting, "--", sim,"/minrowcor_bench.rds"))
  
  bestk_bic_main <- which.min(BICscores_main[BICscores_main!=0])
  bestk_cor_main <- which.min(minrowcor_main[minrowcor_main!=0])
  # bestk_bic_bench <- which.min(BICscores_bench[BICscores_bench!=0])
  # bestk_cor_bench <- which.min(minrowcor_bench[minrowcor_bench!=0])
  for(j in 1:length(threshes)){
    allShdS[[j]] <- get_shd_for_several_methods(kmainbic = bestk_bic_main,
                                           # kbenchbic = bestk_bic_bench,
                                           kmaincor = bestk_cor_main,
                                           # kbenchcor = bestk_cor_bench,
                                           today_date = today_date,
                                           tt = sim,
                                           gesshd = NULL, 
                                           pcshd = NULL,
                                           adjmat_trueCPDAG = adjmat_trueCPDAG,
                                           thresh = threshes[j],
                                           estimands = estimands,
                                           ordered = F)
    cat("[INFO]: j =", j, "\n")
  }  
  saveRDS(allShdS, file = "allshd2.rds")
  
}
