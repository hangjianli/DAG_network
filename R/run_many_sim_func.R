# simulation function -------------------------------------------------------
run_many_sim <- function(args, max.iter = 30, seed = 10, lambda2 = 0.2,
                         lambda.len = 10, threshold = seq(0,0.2,length.out = 4)){
  today_date <- Sys.Date()
# Generate estimands ------------------------------------------------------
  ## If estimands are given, comment out the line below
  # estimands <- generate_parameters(args, seed = seed)
  png("B_true.png")
  # 2. Create the plot
  heatmap.2(abs(estimands$b), dendrogram = "none", Rowv = F, Colv = F, trace = "none")
  # 3. Close the file
  dev.off()
  saveRDS(estimands, file = paste0(format(today_date, "%Y-%m-%d"), "-vers-", args$setting, "-estimands",".rds"))
  saveRDS(args, file = paste0(format(today_date, "%Y-%m-%d"), "-vers-", args$setting, "-args",".rds"))
  p = estimands$realp
  summary <- vector("list", length = args$num_sim)
  for(w in 1:args$num_sim){
    sim_index = paste0(today_date, '-', args$setting, "-", w)
    # # If X is given, comment out the section below
    # X_ <- sim_X(vers = w, p = p,
    #             args = args, omg.sq = estimands$omg.sq,
    #             sig = estimands$sig, b = estimands$b)
    # X <- X_$X
# Run all models once ---------------------------------------------------
    results <- run_three_model(X = X, 
                               vers = w,
                               estimands = estimands, 
                               lambda2 = lambda2,
                               max.iter = max.iter,
                               lambda.len = lambda.len,
                               threshold = threshold)
    saveRDS(results, file = paste0(format(today_date, "%Y-%m-%d"),
                                   "-sim-", w, "-results", ".rds"))
# Gather summary stats ----------------------------------------------------
    summary[[w]] <- calculate_summary_stats(results, estimands, args)
    cat("This is the ", w, "th sample!\n")
  }
  return(list(summary=summary,
              estimands=estimands))
}
