num_sim <- args$num_sim
threshes <- seq(0, 0.8,length.out = 20)
allShdS <- vector(mode = "list", length = length(threshes))

today_date <- Sys.Date()
# today_date <- as.Date("2019-05-16")

for(i in 1:10){
  setwd(paste0("~/Dropbox/research/code/", args$setting,"--", i))
  estimands <- readRDS(paste0("~/Dropbox/research/code/", args$setting,"--", i, "/estimands.rds"))
  

  # X_ <- readRDS(file = paste0(format(today_date, "%Y-%m-%d"), "-vers-",
  #                           args$setting, "-X-",i,".rds"))
  # X <- X_$X
  # n <- dim(X)[1]
  # p <- dim(X)[2]
  # dimnames(estimands$b) <- list(as.character(1:p), as.character(1:p))
  # adjmat_trueCPDAG <- bnstruct::dag.to.cpdag(1*(estimands$b != 0)) 
  # saveRDS(adjmat_trueCPDAG, "adjmat_trueCPDAG.rds")
  
  allshd <- readRDS(paste0("~/Dropbox/research/code/", args$setting,"--", i, "/allshd2.rds"))
  BICscores_bench <- readRDS(paste0("~/Dropbox/research/code/", args$setting,"--", i, "/BICscores_bench.rds"))
  minrowcor_bench <- readRDS(paste0("~/Dropbox/research/code/", args$setting,"--", i, "/minrowcor_bench.rds"))
  BICscores_main <- readRDS(paste0("~/Dropbox/research/code/", args$setting,"--", i, "/BICscores_main.rds"))
  minrowcor_main <- readRDS(paste0("~/Dropbox/research/code/", args$setting,"--", i, "/minrowcor_main.rds"))
  adjmat_trueCPDAG <- readRDS(paste0("~/Dropbox/research/code/", args$setting,"--", i, "/adjmat_true_dag.rds"))
  # adjmat_trueCPDAG <- readRDS(paste0("~/Dropbox/research/code/", args$setting,"--", i, "/adjmat_trueCPDAG.rds"))
  
  bestk_bic_main <- which.min(BICscores_main[BICscores_main!=0])
  bestk_cor_main <- which.min(minrowcor_main[minrowcor_main!=0])
  bestk_bic_bench <- which.min(BICscores_bench[BICscores_bench!=0])
  bestk_cor_bench <- which.min(minrowcor_bench[minrowcor_bench!=0])
  
  # gesshd <- readRDS(paste0("~/Dropbox/research/code/", args$setting,"--", i, "/X-", i,"-shdX_ges.rds"))
  # pcshd <- readRDS(paste0("~/Dropbox/research/code/", args$setting,"--", i, "/shd_pc.rds"))
  # shdL_pc <- readRDS(paste0("~/Dropbox/research/code/", args$setting,"--", i, "/shdL_pc.rds"))
  # shdXLBIC <- readRDS(paste0("~/Dropbox/research/code/", args$setting,"--", i, "/shdXLBIC.rds"))
  
  for(j in 1:length(threshes)){
    allShdS[[j]] <- get_shd_for_several_methods(
                                                kmainbic = bestk_bic_main,
                                                kbenchbic = bestk_bic_bench,
                                                kmaincor = bestk_cor_main,
                                                kbenchcor = bestk_cor_bench,
                                                today_date = today_date,
                                                tt = i,
                                                gesshd = NULL,
                                                pcshd = NULL,
                                                adjmat_trueCPDAG = adjmat_trueCPDAG,
                                                thresh = threshes[j],
                                                # thresh = 0.05,
                                                estimands = estimands,
                                                ordered = T)
    
  # allShdS$shdXLBIC <- unlist(shdXLBIC)
  # allShdS$shdL_pc <- unlist(shdL_pc)   
  # saveRDS(data.frame(allShdS), file = "allshd.rds")
    # allShdS[[i]]$shdXLBIC <- unlist(shdXLBIC)
    # allShdS[[i]]$shdL_pc <- unlist(shdL_pc) 
    cat("[INFO]: j =", j, "\n")
  }
  
  bestbench <- which.min(sapply(allShdS, function(x) abs(x$shdXbench['pnum'] - allshd$shdXmain['pnum'])))
  allshd$shdXbench <- allShdS[[bestbench]]$shdXbench
  allshd$shdXbenchCor <- allShdS[[bestbench]]$shdXbenchCor
  saveRDS(allshd, "allshd2.rds")

  cat("[INFO] Sim ", i, "is done. \n")
  
}

# data.frame(allShdS[[3]])
# ordered -----------------------------------------------------------------

# args <- readRDS("~/OneDrive/sims/4291/args.rds")

# estimands <- readRDS("~/OneDrive/sims/4245AR/estimands.rds")
# allShdS <- readRDS("~/OneDrive/sims/4245AR/allShdS.rds")
# setwd(paste0("~/Dropbox/research/code/", args$setting))
# saveRDS(allShdS, "allShdS.rds")

# png(paste0(args$setting, "-4.png"),width = 500, height = 460)
    


setwd(paste0("~/Dropbox/research/code/", args$setting))
p <- estimands$realp
n <- args$n
rocmain <- sapply(allShdS, function(x) c(FP =(x$shdXmain['pnum'] - x$shdXmain['TP']),
                                         TP = x$shdXmain['TP'] ))
rocbench <- sapply(allShdS, function(x) c(FP =(x$shdXbench['pnum'] - x$shdXbench['TP']),
                                          TP = x$shdXbench['TP']))
xmax <- max(rocmain[1,])*7/8
ymax <- max(rocmain[2,])

setEPS()
postscript(paste0(args$setting, ".eps"))
mar.default <- c(4,4,2,2)
par(mar = mar.default + c(1, 1, 0, 0)) 
plot(t(rocmain),
     xlim=c(0,xmax),
     ylim = c(50, ymax),
     type = "l",
     lty = 4,
     lwd = 4,
     col="red",
     xlab = "False Positive",
     ylab= "True Positive",
     xaxt = "n",
     yaxt = "n",
     cex.axis = 2,
     cex.lab = 2)
points(t(rocbench),
       col = "blue",
       lty = 5,
       lwd = 4,
       type="l")
axis(side = 1, at = seq(0, xmax, length.out = 5),
     labels = round(seq(0, xmax, length.out = 5),1),
     cex.axis = 1.5)
axis(side = 2, at = seq(50, ymax,length.out = 5), 
     labels = seq(50, ymax,length.out = 5), cex.axis = 1.5)
legend("bottomright", legend = c("BCD", "bench"),
       col = c("red", "blue"), lty = c(4,5), lwd = c(4,4),
       cex = 2)
# legend("bottomright", legend = c("BCD"),
#        col = c("red"), lty = c(4))
dev.off()

