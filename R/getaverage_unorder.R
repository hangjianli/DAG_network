options(digits = 3)
nsim = as.numeric(args$num_sim)
# nsim = 5
allshd <- readRDS(paste0("~/Dropbox/research/code/", as.character(args$setting), '--1/', "allshd.rds"))
setting <- as.character(args$setting)
num_statistic <- length(rownames(allshd)) # number of summary statistics (rows in table)
total <- data.frame(row.names = rownames(allshd), 
                    shdXmain=rep(0, num_statistic), 
                    shdXmainCor=rep(0, num_statistic),
                    shdGES=rep(0, num_statistic),
                    shdPC = rep(0, num_statistic),
                    GESL = rep(0, num_statistic),
                    PCL = rep(0, num_statistic),
                    sb = rep(0, num_statistic),
                    sbL = rep(0, num_statistic))



theta_re <- rep(0, nsim)
Bre_average <- numeric(length = 5)
thetaRE_average <- numeric(length = 3)
for (i in 1:args$num_sim){
  #compute relative error
  # Bre <- readRDS(paste0("~/Dropbox/research/code/",setting, "--", i, "/Bre.rds"))
  # Bre_average <- Bre_average + unlist(Bre)
  Thetare <- readRDS(paste0("~/Dropbox/research/code/",setting, "--", i, "/thetaRE.rds"))
  thetaRE_average <- thetaRE_average + unlist(Thetare)
  allshd <- readRDS(paste0("~/Dropbox/research/code/",setting, "--", i, "/allshd.rds"))  
  # allshd <- as.data.frame(allshd)
  # allshd$shdXflipbic <- unlist(shdXflipbic)
  # allshd <- allshd[,names(total)]
  total <- total + allshd
  # theta_re[i] <- readRDS(paste0("~/Dropbox/research/code/",setting, "--", i, "/thetaRE.rds"))
}

# saveRDS(theta_re, paste0("../", setting,"/thetaER.rds"))
Bre_average <- Bre_average / nsim
thetaRE_average <- thetaRE_average / nsim

total <- round(total / nsim, 2)
total$B0 = c(estimands$s0, rep(0,num_statistic-1))
saveRDS(total, paste0("../", setting,"/total.rds"))
saveRDS(Bre_average, paste0("../", setting,"/BreTotal.rds"))
saveRDS(thetaRE_average, file = paste0("../", setting,"/thetaRETotal.rds"))


# boxplot -----------------------------------------------------------------


# shdmain <- shdges <-  shdpc <- shdgesL <- shdpcL <- rep(0, 30)
shdgesdiff <- shdpcdiff <- shdsbdiff <- rep(0,4*10)
JIgesdiff <- JIpcdiff <- JIsbdiff <- rep(0,4*10)
  
# shddiff <- rep(0,10)
# shdgesdiff <- shdgesdiff[1:40]
# shdpcdiff <- shdpcdiff[1:40]
  
args <- readRDS("~/Dropbox/research/code/072630/args.rds")
setting <- as.character(args$setting)
# settting <- "4305"
nsim <- 10
a <- 30
for (i in 1:nsim){
  #compute relative error
  # Bre <- readRDS(paste0("~/Dropbox/research/code/",setting, "--", i, "/Bre.rds"))
  # Bre_average <- Bre_average + unlist(Bre)
  # Thetare <- readRDS(paste0("~/Dropbox/research/code/",setting, "--", i, "/thetaRE.rds"))
  # thetaRE_average <- thetaRE_average + unlist(Thetare)
  allshd <- readRDS(paste0("~/Dropbox/research/code/", setting, "--", i, "/allshd.rds"))  
  # allshd <- as.data.frame(allshd)
  # shdmain[as.integer(i+a)] <- allshd$shdXmain[7]
  # shdges[as.integer(i+a)] <- allshd$shdGES[7]
  # shdpc[as.integer(i+a)] <- allshd$shdPC[7]
  # shdpcL[as.integer(i+a)] <- allshd$shdL_pc[7]
  # shdgesL[as.integer(i+a)] <- allshd$shdXLBIC[7]
  # print(c(allshd$shdGES[7], allshd$shdXLBIC[7]))
  shdgesdiff[as.integer(i+a)] <- allshd$shdGES[7] - allshd$shdXLBIC[7]
  shdpcdiff[as.integer(i+a)] <- allshd$shdPC[7] - allshd$shdL_pc[7]
  shdsbdiff[as.integer(i+a)] <- allshd$shdSB[7] - allshd$shdSBL[7]
  
  JIgesdiff[as.integer(i+a)] <- allshd$shdXLBIC[5] - allshd$shdGES[5]
  JIpcdiff[as.integer(i+a)] <- allshd$shdL_pc[5] - allshd$shdPC[5]
  JIsbdiff[as.integer(i+a)] <- allshd$shdSBL[5] - allshd$shdSB[5]
  
  
  # shddiff[as.integer(i+a)] <- allshd$shdGES[6] - allshd$shdXLBIC[6]
  
  # GESL <- readRDS(paste0("~/Dropbox/research/code/", setting, "--", i, "/shdXLBIC.rds"))
  # PCL <- readRDS(paste0("~/Dropbox/research/code/", setting, "--", i, "/shdL_pc.rds"))
  # allshd$GESL <- unlist(GESL)
  # allshd$PCL <- PCL
  # shdXflipbic <- readRDS(paste0("~/OneDrive/Documents/research/code/11072--", i, "/shdXX_flip_bic.rds"))
  # allshd$shdXflipbic <- unlist(shdXflipbic)
  # theta_re[i] <- readRDS(paste0("~/Dropbox/research/code/",setting, "--", i, "/thetaRE.rds"))
}
# boxplot(shddiff)
# mydata <- data.frame(shd = c(shdmain, shdges, shdgesL, shdpc, shdpcL),
#                      id = rep(rep(c("andes", "hailfinder", "barley"), each=args$num_sim), 5),
#                      label = rep(c("BCD","GES", "GESL", "PC", "PCL"), each= 30))

# mydata_ordered <- data.frame(shd = shddiff,
#                              id = rep(rep(c("arth150", "magic-irri", "ecoli70"), each=args$num_sim), 2),
#                              label = rep(c("facebook", "celegans_n306"), each= 30))

mydata <- data.frame(shddiff = c(shdgesdiff, shdpcdiff, shdsbdiff),
                     id = rep(rep(rep(c("arth150", "ecoli70"), each=args$num_sim), 2),3),
                     label = rep(c("GES", "PC", "SPB"), each = 40),
                     net = rep(rep(c("facebook", "celegans_n306"), each = 20),3))


mydataSHD <- data.frame(shddiff = c(shdgesdiff, shdpcdiff),
                     id = rep(rep(rep(c("Toeplitz", "equi.cor", "AR", "star"), each=args$num_sim), 2), 2),
                     label = rep(rep(c("n<p", "n>p"), each = 40),2),
                     method = rep(c("GES", "PC"), each = 80))


# saveRDS(mydata, "unordered_sim_SHD_JI.rds")

dataToep <- subset(mydataSHD, id == "Toeplitz")
dataEqui <- subset(mydataSHD, id == "equi.cor")
dataAR <- subset(mydataSHD, id == "AR")
datastar <- subset(mydataSHD, id == "star")



# png(paste0("AR", ".png"), width = 500, height = 460)
setEPS()
postscript(paste0("Toeplitz", ".eps"))
par(mar = c(2,2,2,2))
ggplot(dataToep, aes(x = label, y = shddiff)) + 
  geom_boxplot(aes(fill = method), width=.3) + 
  xlab("Sample size (n, p)") + 
  ylab("Decrease in SHD") + 
  # facet_grid(.~label) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 1.5) + 
  theme(
    strip.text.x = element_text(size = 18),
    panel.background = element_blank(),
    panel.grid.minor = element_blank(),
    # panel.border = element_rect(color = "black", fill = NA),
    axis.line = element_line(colour = "black"),
    axis.title=element_text(size=25,face="bold"),
    axis.text = element_text(size=25),
    legend.text=element_text(size=20),
    # axis.text.x=element_blank(),
    legend.title = element_text(size = 20),
    strip.background = element_rect(color = "black", fill = "white")
  ) + 
  scale_x_discrete(labels=c("n<p" = "(100, 200)", "n>p" = "(300, 100)"))
dev.off()






# Jaccard -----------------------------------------------------------------



mydataJI <- data.frame(JIdiff = c(JIgesdiff, JIpcdiff, JIsbdiff),
                       id = rep(rep(rep(c("arth150", "ecoli70"), each=args$num_sim), 2),3),
                     label = rep(c("GES", "PC", "SPB"), each = 40),
                     net = rep(rep(c("facebook", "celegans_n306"), each = 20),3))

saveRDS(mydataJI, "JIdiff_last.rds")



# png(paste0("SHDdiff", ".png"), width = 500, height = 460)
setEPS()
postscript("JIdiffnew.eps")
par(mar = c(2,2,2,2))
gplot <-
  ggplot(mydataJI, aes(x = id, y = JIdiff)) + 
  geom_boxplot(aes(fill = label)) + 
  xlab("Real DAG networks") + 
  ylab("Increase in Jaccard index") + 
  facet_grid(.~net) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 1.5) + 
  theme(
    strip.text.x = element_text(size = 48),
    panel.background = element_blank(),
    panel.grid.minor = element_blank(),
    # panel.border = element_rect(color = "black", fill = NA),
    axis.line = element_line(colour = "black"),
    axis.title=element_text(size=48,face="bold"),
    axis.text = element_text(size=48),
    legend.text=element_text(size=48),
    legend.title = element_text(size = 33),
    strip.background = element_rect(color = "black", fill = "white")
  ) +
  scale_y_continuous(limits = c(-0.02,0.13))
ggsave(gplot, filename = "JIdiff_new.eps", device = "eps")
dev.off()

gplot <- ggplot(mydata, aes(x = id, y = shddiff)) + 
  geom_boxplot(aes(fill = label)) + 
  xlab("Real DAG networks") + 
  ylab("Decrease in SHD") + 
  facet_grid(.~net) +
  scale_y_continuous(limits = c(-40,100)) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 1.5) + 
  theme(
    strip.text.x = element_text(size = 48),
    panel.background = element_blank(),
    panel.grid.minor = element_blank(),
    # panel.border = element_rect(color = "black", fill = NA),
    axis.line = element_line(colour = "black"),
    axis.title=element_text(size=48,face="bold"),
    axis.text = element_text(size=48),
    legend.text=element_text(size=48),
    legend.title = element_text(size = 33),
    strip.background = element_rect(color = "black", fill = "white")
  ) 
  # scale_y_continuous(limits = c(-50,100))
ggsave(gplot, filename = "SHDdiff.eps", device = "eps")






# don’t touch -------------------------------------------------------------
mydata_facebook_unordered <- mydata
# saveRDS(mydata_facebook_unordered, "mydata_facebook_unordered.rds")
mydata_celegans_n306_unordered <- mydata
# saveRDS(mydata_celegans_n306_unordered, "mydata_celegans_n306_unordered.rds")
