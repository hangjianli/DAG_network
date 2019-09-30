calculate_summary_stats <- function(results, estimands, args){
  threshold <- results$threshold
  lambda.best.temp <- vector("list", length = length(threshold))
  bench <- main <- twostep <- 
    bench.star <- twostep.star <- 
    data.frame(Threshold = threshold, 
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
  # For each threshold value, calculate best lambda w.r.t. hamming d --------------------------------
  for(i in 1:length(threshold)) {
    SHD <- data.frame(
      twostep = sapply(results$twoStep.hamming, "[[", i),
      bench = sapply(results$bench.hamming, "[[", i),
      main = sapply(results$main.hamming, "[[", i))
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
    bench$lambda.pos[i] <- round(lambda.best.temp[[i]]["bench"],0)
    bench$P[i] <- results$bench.edgenum[i, bench$lambda.pos[i]]
    bench$FP[i] <- results$bench_TPFP[[bench$lambda.pos[i]]]$FP[i]
    bench$TP[i] <- results$bench_TPFP[[bench$lambda.pos[i]]]$TP[i]
    bench$FDR[i] <- results$bench_TPFP[[bench$lambda.pos[i]]]$FP[i]/(results$bench_TPFP[[bench$lambda.pos[i]]]$FP[i] + 
                                                                            results$bench_TPFP[[bench$lambda.pos[i]]]$TP[i])
    bench$SHD[i] <- results$bench.hamming[[bench$lambda.pos[i]]][i]
    bench$JI[i] <- bench$TP[i] / (bench$P[i] + sum(abs(estimands$b) > 1e-6) - bench$TP[i])
    bench$Bl2[i] <- results$berror_bench[i, bench$lambda.pos[i]]
    bench$TE[i] <- results$theta_errors[bench$lambda.pos[i], "lassoIdent"]
    #main -------------------
    main$lambda.pos[i] <- round(lambda.best.temp[[i]]["main"],0)
    main$P[i] <- results$main.edgenum[i, main$lambda.pos[i]]
    main$FP[i] <- results$main_TPFP[[main$lambda.pos[i]]]$FP[i]
    main$TP[i] <- results$main_TPFP[[main$lambda.pos[i]]]$TP[i]
    main$FDR[i] <- results$main_TPFP[[main$lambda.pos[i]]]$FP[i]/(results$main_TPFP[[main$lambda.pos[i]]]$FP[i] + 
                                                                      results$main_TPFP[[main$lambda.pos[i]]]$TP[i])
    main$SHD[i] <- results$main.hamming[[main$lambda.pos[i]]][i]
    main$JI[i] <- main$TP[i] / (main$P[i] + sum(abs(estimands$b) > 1e-6) - main$TP[i])
    main$TE[i] <- results$theta_errors[main$lambda.pos[i], "main"]
    main$Bl2[i] <- results$berror_main[i, main$lambda.pos[i]]
    #twostep -------------------
    twostep$lambda.pos[i] <- round(lambda.best.temp[[i]]["twostep"], 0)
    twostep$P[i] <- results$twostep.edgenum[i, twostep$lambda.pos[i]]
    twostep$FP[i] <- results$twoStep_TPFP[[twostep$lambda.pos[i]]]$FP[i]
    twostep$TP[i] <- results$twoStep_TPFP[[twostep$lambda.pos[i]]]$TP[i]
    twostep$FDR[i] <- results$twoStep_TPFP[[twostep$lambda.pos[i]]]$FP[i]/(results$twoStep_TPFP[[twostep$lambda.pos[i]]]$FP[i] + 
                                                                       results$twoStep_TPFP[[twostep$lambda.pos[i]]]$TP[i])
    twostep$SHD[i] <- results$twoStep.hamming[[twostep$lambda.pos[i]]][i]
    twostep$JI[i] <- twostep$TP[i] / (twostep$P[i] + sum(abs(estimands$b) > 1e-6) - twostep$TP[i])
    twostep$TE[i] <- results$theta_errors[twostep$lambda.pos[i], "scaledLasso"]
    twostep$Bl2[i] <- results$berror_twoStep[i, twostep$lambda.pos[i]]
  }
  
  output[seq(1,by = 5,length.out = length(threshold)), 3:11] <- main[, -1]
  output[seq(2, by = 5, length.out = length(threshold)), 3:11] <- bench[,-1]
  output[seq(4, by = 5, length.out = length(threshold)), 3:11] <- twostep[,-1]
  
  bench_num <- mapply(function(x, y) which.min(abs(x-y)), bench$P, as.data.frame(t(results$main.edgenum)))
  twostep_num <- mapply(function(x, y) which.min(abs(x-y)), twostep$P, as.data.frame(t(results$main.edgenum)))
  bench.star$lambda.pos <- bench_num
  twostep.star$lambda.pos <- twostep_num
  
  for(i in 1:length(threshold)){
    # match for bench --------------------------------------------------------
    bench.star$P[i] <- results$main.edgenum[i, bench.star$lambda.pos[i]]
    bench.star$FP[i] <- results$main_TPFP[[bench.star$lambda.pos[i]]]$FP[i]
    bench.star$TP[i] <- results$main_TPFP[[bench.star$lambda.pos[i]]]$TP[i]
    bench.star$FDR[i] <- results$main_TPFP[[bench.star$lambda.pos[i]]]$FP[i]/(results$main_TPFP[[bench.star$lambda.pos[i]]]$FP[i] + 
                                                                          results$main_TPFP[[bench.star$lambda.pos[i]]]$TP[i])
    bench.star$SHD[i] <- results$main.hamming[[bench.star$lambda.pos[i]]][i]
    bench.star$JI[i] <- bench.star$TP[i] / (bench.star$P[i] + sum(abs(estimands$b) > 1e-6) - bench.star$TP[i])
    bench.star$TE[i] <- results$theta_errors[bench.star$lambda.pos[i], "main"]
    bench.star$Bl2[i] <- results$berror_main[i, bench.star$lambda.pos[i]]
    # Match for twostep -------------------------------------------------------
    twostep.star$P[i] <- results$main.edgenum[i, twostep.star$lambda.pos[i]]
    twostep.star$FP[i] <- results$main_TPFP[[twostep.star$lambda.pos[i]]]$FP[i]
    twostep.star$TP[i] <- results$main_TPFP[[twostep.star$lambda.pos[i]]]$TP[i]
    twostep.star$FDR[i] <-  results$main_TPFP[[twostep.star$lambda.pos[i]]]$FP[i]/(results$main_TPFP[[twostep.star$lambda.pos[i]]]$FP[i] + 
                                                                                   results$main_TPFP[[twostep.star$lambda.pos[i]]]$TP[i])
    twostep.star$SHD[i] <- results$main.hamming[[twostep.star$lambda.pos[i]]][i]
    twostep.star$JI[i] <-  twostep.star$TP[i] / (twostep.star$P[i] + sum(abs(estimands$b) > 1e-6) - twostep.star$TP[i])
    twostep.star$TE[i] <- results$theta_errors[twostep.star$lambda.pos[i], "main"]
    twostep.star$Bl2[i] <- results$berror_main[i, twostep.star$lambda.pos[i]]
  }
  
  output[seq(3, by = 5, length.out = length(threshold)), 3:11] <- bench.star[,-1]
  output[seq(5, by = 5, length.out = length(threshold)), 3:11] <- twostep.star[,-1]
  
  return(output)
}
