options(digits = 3)
nsim = as.numeric(args$num_sim)
# nsim <- 5

allshd <- readRDS(paste0("~/Dropbox/research/code/", as.character(args$setting), '--1/', "allshd2.rds"))
setting <- as.character(args$setting)
num_statistic <- length(names(allshd[[1]])) # number of summary statistics (rows in table)
total <- data.frame(row.names = names(allshd[[1]]), 
                    shdXmain=rep(0, num_statistic),
                    shdXbench=rep(0, num_statistic),
                    shdXmainCor=rep(0, num_statistic),
                    shdXbenchCor=rep(0, num_statistic)
                    # , 
                    
                    # , 
                    # shdGES=rep(0, num_statistic),
                    # shdPC = rep(0, num_statistic), 
                    # GESL = rep(0, num_statistic),
                    # PCL = rep(0, num_statistic)
                    # shdXflip = rep(0,num_statistic)
                    # shdXflipcor = rep(0,num_statistic)
                    )



theta_re <- rep(0, nsim)
Bre_average <- numeric(length = 5)
thetaRE_average <- numeric(length = 3)
for (i in 1:nsim){
  #compute relative error
  # Bre <- readRDS(paste0("~/Dropbox/research/code/",setting, "--", i, "/Bre.rds"))
  # Bre_average <- Bre_average + unlist(Bre)
  # Thetare <- readRDS(paste0("~/Dropbox/research/code/",setting, "--", i, "/thetaRE.rds"))
  # thetaRE_average <- thetaRE_average + unlist(Thetare)
  allshd <- readRDS(paste0("~/Dropbox/research/code/",setting, "--", i, "/allshd2.rds"))  
  # allshd <- as.data.frame(allshd)
  allshd <- as.data.frame(allshd[1:4])
  # GESL <- readRDS(paste0("~/Dropbox/research/code/", setting, "--", i, "/shdXLBIC.rds"))
  # PCL <- readRDS(paste0("~/Dropbox/research/code/", setting, "--", i, "/shdL_pc.rds"))
  # allshd$GESL <- unlist(GESL)
  # allshd$PCL <- PCL
  # shdXflipbic <- readRDS(paste0("~/OneDrive/Documents/research/code/11072--", i, "/shdXX_flip_bic.rds"))
  # allshd$shdXflipbic <- unlist(shdXflipbic)
  allshd <- allshd[,names(total)]
  total <- total + allshd
  # theta_re[i] <- readRDS(paste0("~/Dropbox/research/code/",setting, "--", i, "/thetaRE.rds"))
}

# saveRDS(theta_re, paste0("../", setting,"/thetaER.rds"))
Bre_average <- Bre_average / nsim
thetaRE_average <- thetaRE_average / nsim

total <- round(total / nsim, 2)
total$B0 = c(estimands$s0, rep(0,num_statistic-1))
saveRDS(total, paste0("../", setting,"/total2.rds"))
saveRDS(Bre_average, paste0("../", setting,"/BreTotal.rds"))
saveRDS(thetaRE_average, file = paste0("../", setting,"/thetaRETotal.rds"))



# boxplot -----------------------------------------------------------------

shddiff <- JIdiff <-  rep(0, 60)
shddiff <- rep(0,10)

# args <- readRDS("~/Dropbox/research/code/4305/args.rds")
setting <- as.character(args$setting)
# settting <- "4305"
nsim <- 10
a <- 50
for (i in 1:nsim){
  #compute relative error
  # Bre <- readRDS(paste0("~/Dropbox/research/code/",setting, "--", i, "/Bre.rds"))
  # Bre_average <- Bre_average + unlist(Bre)
  # Thetare <- readRDS(paste0("~/Dropbox/research/code/",setting, "--", i, "/thetaRE.rds"))
  # thetaRE_average <- thetaRE_average + unlist(Thetare)
  allshd <- readRDS(paste0("~/Dropbox/research/code/",setting, "--", i, "/allshd2.rds"))  
  # allshd <- as.data.frame(allshd)
  # shdmain[as.integer(i+a)] <- allshd$shdXmain['myshd']
  # shdbench[as.integer(i+a)] <- allshd$shdXbench['myshd']
  shddiff[as.integer(i+a)] <- allshd$shdXbench['myshd'] - allshd$shdXmain['myshd']
  JIdiff[as.integer(i+a)] <- allshd$shdXmain['JI'] - allshd$shdXbench['JI']
  # GESL <- readRDS(paste0("~/Dropbox/research/code/", setting, "--", i, "/shdXLBIC.rds"))
  # PCL <- readRDS(paste0("~/Dropbox/research/code/", setting, "--", i, "/shdL_pc.rds"))
  # allshd$GESL <- unlist(GESL)
  # allshd$PCL <- PCL
  # shdXflipbic <- readRDS(paste0("~/OneDrive/Documents/research/code/11072--", i, "/shdXX_flip_bic.rds"))
  # allshd$shdXflipbic <- unlist(shdXflipbic)
  # theta_re[i] <- readRDS(paste0("~/Dropbox/research/code/",setting, "--", i, "/thetaRE.rds"))
}

# mydata <- data.frame(shd = c(shdmain, shdbench),
#                      id = rep(c("andes", "hailfinder", "barley"), each=args$num_sim),
#                      label = rep(c("BCD", "Bench"), each= 30))



real_data_ordered <- data.frame(shddiff = shddiff,
                                id = rep(rep(c("andes", "hailfinder", "barley"), each=args$num_sim),2),
                                label = rep(c("facebook", "celegans_n306"), each = 30))

ggplot(real_data_ordered, aes(x = id, y = shddiff)) + 
  geom_boxplot(aes(fill = label)) + 
  xlab("Real DAG networks") + 
  ylab("Decrease in SHD") + 
  # facet_grid(.~net) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 1.5) + 
  theme(
    strip.text.x = element_text(size = 18),
    panel.background = element_blank(),
    panel.grid.minor = element_blank(),
    # panel.border = element_rect(color = "black", fill = NA),
    axis.line = element_line(colour = "black"),
    axis.title=element_text(size=17,face="bold"),
    axis.text = element_text(size=16),
    legend.text=element_text(size=13),
    legend.title = element_text(size = 13),
    strip.background = element_rect(color = "black", fill = "white")
  )

png(paste0("celegans_n306_ordered", ".png"))
ggplot(mydata, aes(x = id, y = shd)) + 
  geom_boxplot(aes(fill = label)) +
  xlab("Real DAG networks") + 
  ylab("SHD")
dev.off()






# don’t touch -------------------------------------------------------------
mydata_facebook_ordered <- mydata
saveRDS(mydata_facebook_ordered, "mydata_facebook_ordered.rds")
mydata_celegans_n306 <- mydata
