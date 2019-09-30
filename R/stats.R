####################################
#####
#####  
#####
#####
##### input: a table with best lambdas for each threshold "lambda.best.temp"
##### output: key stats in data frame
#####
########################################################################


################################
###### some initializations ####
################################

lambda.best.temp <- vector("list", length = length(threshold))

bench <- main <- twostep <- bench.star <- twostep.star <- data.frame(Threshold = threshold, 
                                                                     lambda.pos = rep(0,length(threshold)),
                                                                     P = rep(0,length(threshold)),
                                                                     FP = rep(0,length(threshold)), 
                                                                     TP = rep(0,length(threshold)),  
                                                                     FDR = rep(0,length(threshold)),  
                                                                     SHD = rep(0,length(threshold)),
                                                                     JI = rep(0,length(threshold)),
                                                                     TE = rep(NaN,length(threshold)),
                                                                     Bl2 = rep(0,length(threshold)))

output <- data.frame(Threshold = unlist(lapply(threshold, rep, 5)),
                     Method = rep(c("Main", "Bench", "Bench*", "Twostep", "Twostep*"), length(threshold)),
                     lambda.pos = rep(0, 5*length(threshold)),
                     P = rep(0, 5* length(threshold)),
                     FP = rep(0, 5* length(threshold)), 
                     TP = rep(0, 5*length(threshold)),  FDR = rep(0, 5*length(threshold)),  
                     SHD = rep(0, 5*length(threshold)),
                     JI = rep(0, 5*length(threshold)), TE = rep(0, 5*length(threshold)),
                     BE = rep(0, 5*length(threshold))
                     ) 

#####################################################################################
##### For each threshold value, calculate best lambda w.r.t. hamming distance #######
#####################################################################################

for(i in 1:length(threshold)) {
  SHD <- data.frame(
    twostep = sapply(twoStep.hamming, "[[", i),
    bench = sapply(bench.hamming, "[[", i),
    main = sapply(main.hamming, "[[", i))
  lambda.best.temp[[i]] <- sapply(SHD, which.min)
}

#####################################################################################
####### For each threshold value, calculate stats for the best lambda ##############
####### P: number of estimated edges
####### FP: false positive edges
####### TP: true positive  edges
####### FDR: false discovery rate
####### SHD: structural hamming distance
####### JI: jaccard index
#####################################################################################

 
for(i in 1:length(threshold)){
  #benchmark---------------------
  # bbench <- readRDS(paste0(n,p, "_", lambda.best.temp[[i]]["bench"], "_b_bench.rds"))
  bbench <- readRDS(paste0(n,p, "_", lambda.best.temp[[i]]["bench"], "_b_bench.rds"))
  bbench <- ifelse(abs(bbench) > threshold[i], sign(bbench)*(abs(bbench) - threshold[i]), 0)
  bench$lambda.pos[i] <- round(lambda.best.temp[[i]]["bench"],0)
  bench$P[i] <- sum(abs(bbench) > 1e-7)
  bench$FP[i] <- bench_TPFP[[lambda.best.temp[[i]]["bench"]]][[1]][i]
  bench$TP[i] <-   bench_TPFP[[lambda.best.temp[[i]]["bench"]]][[2]][i]
  bench$FDR[i] <- bench_TPFP[[lambda.best.temp[[i]]["bench"]]][[1]][i]/(bench_TPFP[[lambda.best.temp[[i]]["bench"]]][[1]][i] + 
                                                                          bench_TPFP[[lambda.best.temp[[i]]["bench"]]][[2]][i])
  bench$SHD[i] <- bench.hamming[[lambda.best.temp[[i]]["bench"]]][[i]]
  bench$JI[i] <-  bench$TP[i] / (bench$P[i] + sum(abs(b) > 1e-6) - bench$TP[i])
  bench$Bl2[i] <- norm(bbench - b, "f")/s0
  bench$TE[i] <- theta_errors[bench$lambda.pos[i], "lassoIdent"]
  #main -------------------
  bmain <- readRDS(paste0(n,p, "_", lambda.best.temp[[i]]["main"], "_b_main.rds"))
  bmain <- ifelse(abs(bmain) > threshold[i], sign(bmain)*(abs(bmain) - threshold[i]), 0)
  main$lambda.pos[i] <- round(lambda.best.temp[[i]]["main"],0)
  main$P[i] <- sum(abs(bmain) > 1e-7)
  main$FP[i] <-  main_TPFP[[lambda.best.temp[[i]]["main"]]][[1]][i]
  main$TP[i] <- main_TPFP[[lambda.best.temp[[i]]["main"]]][[2]][i]
  main$FDR[i] <- main_TPFP[[lambda.best.temp[[i]]["main"]]][[1]][i]/(main_TPFP[[lambda.best.temp[[i]]["main"]]][[1]][i] + 
                                                                       main_TPFP[[lambda.best.temp[[i]]["main"]]][[2]][i])
  main$SHD[i] <- main.hamming[[lambda.best.temp[[i]]["main"]]][[i]]
  main$JI[i] <- main$TP[i] / (main$P[i] + sum(abs(b) > 1e-6) - main$TP[i])
  main$TE[i] <- theta_errors[main$lambda.pos[i], "main"]
  main$Bl2[i] <- norm(bmain - b, "f") / s0
  #twostep -------------------
  btwostep <- readRDS(paste0(n,p, "_", lambda.best.temp[[i]]["twostep"], "_b_twostep.rds"))
  btwostep <- ifelse(abs(btwostep) > threshold[i], sign(btwostep)*(abs(btwostep) - threshold[i]), 0)
  twostep$lambda.pos[i] <- round(lambda.best.temp[[i]]["twostep"], 0)
  twostep$P[i] <- sum(abs(btwostep) > 1e-7)
  twostep$FP[i] <- twoStep_TPFP[[lambda.best.temp[[i]]["twostep"]]][[1]][i]
  twostep$TP[i] <- twoStep_TPFP[[lambda.best.temp[[i]]["twostep"]]][[2]][i]
  twostep$FDR[i] <- twoStep_TPFP[[lambda.best.temp[[i]]["twostep"]]][[1]][i] / (twoStep_TPFP[[lambda.best.temp[[i]]["twostep"]]][[1]][i] +
                                                                                  twoStep_TPFP[[lambda.best.temp[[i]]["twostep"]]][[2]][i] )
  twostep$SHD[i] <- twoStep.hamming[[lambda.best.temp[[i]]["twostep"]]][[i]]
  twostep$JI[i] <- twostep$TP[i] / (twostep$P[i] + sum(abs(b) > 1e-6) - twostep$TP[i])
  twostep$TE[i] <- theta_errors[twostep$lambda.pos[i],"scaledLasso"]
  twostep$Bl2[i] <- norm(btwostep - b, "f") / s0
}

output[seq(1,by = 5,length.out = length(threshold)), 3:11] <- main[, -1]
output[seq(2, by = 5, length.out = length(threshold)), 3:11] <- bench[,-1]
output[seq(4, by = 5, length.out = length(threshold)), 3:11] <- twostep[,-1]
main_num_edge <- matrix(0, length(threshold), length(lambda.path))

for(i in 1:length(lambda.path)){
  main_btemp <- readRDS(paste0(n,p, "_", i, "_b_main.rds"))
  main_num_edge[, i] <- unlist(lapply(lapply(threshold, function(x) ifelse(abs(main_btemp) > x, sign(main_btemp)*(abs(main_btemp) - x), 0)),
                                      function(x) sum(abs(x) > 1e-7)))
}

bench_num <- mapply(function(x, y) which.min(abs(x-y)), bench$P, as.data.frame(t(main_num_edge)))
twostep_num <- mapply(function(x, y) which.min(abs(x-y)), twostep$P, as.data.frame(t(main_num_edge)))
bench.star$lambda.pos <- bench_num
twostep.star$lambda.pos <- twostep_num

for(i in 1:length(threshold)){
  bench.btemp <- readRDS(paste0(n,p, "_", bench.star$lambda.pos[i], "_b_main.rds"))
  bench.btemp <- ifelse(abs(bench.btemp) > threshold[i], sign(bench.btemp)*(abs(bench.btemp) - threshold[i]), 0)
  bench.star$P[i] <- sum(abs(bench.btemp) > 1e-7)
  bench.star$FP[i] <- main_TPFP[[bench.star$lambda.pos[i]]][[1]][i]
  bench.star$TP[i] <- main_TPFP[[bench.star$lambda.pos[i]]][[2]][i]
  bench.star$FDR[i] <- main_TPFP[[bench.star$lambda.pos[i]]][[1]][i] / (main_TPFP[[bench.star$lambda.pos[i]]][[1]][i] + 
                                                                           main_TPFP[[bench.star$lambda.pos[i]]][[2]][i])
  bench.star$SHD[i] <- main.hamming[[bench.star$lambda.pos[i]]][[i]]
  bench.star$JI[i] <- bench.star$TP[i] / (bench.star$P[i] + sum(abs(b) > 1e-7) - bench.star$TP[i]) 
  bench.star$TE[i] <- theta_errors[bench.star$lambda.pos[i], "main"]
  bench.star$Bl2[i] <- norm(bench.btemp - b, "f") / s0
  
  twostep.btemp <- readRDS(paste0(n,p, "_", twostep.star$lambda.pos[i], "_b_main.rds"))
  twostep.btemp <- ifelse(abs(twostep.btemp) > threshold[i], sign(twostep.btemp)*(abs(twostep.btemp) - threshold[i]), 0)
  twostep.star$P[i] <- sum(abs(twostep.btemp) > 1e-7)
  twostep.star$FP[i] <- main_TPFP[[twostep.star$lambda.pos[i]]][[1]][i]
  twostep.star$TP[i] <- main_TPFP[[twostep.star$lambda.pos[i]]][[2]][i]
  twostep.star$FDR[i] <- main_TPFP[[twostep.star$lambda.pos[i]]][[1]][i] / (main_TPFP[[twostep.star$lambda.pos[i]]][[1]][i] + 
                                                                           main_TPFP[[twostep.star$lambda.pos[i]]][[2]][i])
  twostep.star$SHD[i] <- main.hamming[[twostep.star$lambda.pos[i]]][[i]]
  twostep.star$JI[i] <-  twostep.star$TP[i] / (twostep.star$P[i] + sum(abs(b) > 1e-7) - twostep.star$TP[i])
  twostep.star$TE[i] <- theta_errors[twostep.star$lambda.pos[i],"scaledLasso"]
  twostep.star$Bl2[i] <- norm(twostep.btemp - b, "f") / s0
}

output[seq(3, by = 5, length.out = length(threshold)), 3:11] <- bench.star[,-1]
output[seq(5, by = 5, length.out = length(threshold)), 3:11] <- twostep.star[,-1]
# output
saveRDS(output, paste0("Version", vers, "-",n, "-", p,".rds"))

spec <- data.frame(p = p, 
                   n = n,
                   s0 = s0,
                   # cs = clique_size, bmag = b.magnitude,
                   # proB = prob_B, sparsityTheta = sparse.prob, fixzero = fix.zero,
                   # epsilon = epsilon, thetaType = theta_type,
                   dir = test_num)
write.table(spec, "spec.txt", row.names = F, sep = "\t")
saveRDS(output, file = paste0(test_num, ".rds"))
print(xtable(output, digits=c(0,4,0,0,0,0,0,3,0,3,3,6), latex.environment = center),
      hline.after = c(-1, 0, 5, 10, 15, 20), include.rownames = F)


