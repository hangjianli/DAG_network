process_output_ordered_tmp <- function(
  simID='001', 
  estimands,
  args,
  start=1,
  num_sim=args$num_sim,
  thr=0.1
){
  setwd(paste0('output/', simID))
  # args <- readRDS('args.rds')
  # estimands <- readRDS('estimands.rds')
  bstar_adj <- 1*(abs(estimands$b) > 0)
  for(sim in start:num_sim){
    # BCD results
    BICscores_main <- readRDS(paste0(simID, '--', sim, '/BICscores_main.rds'))  
    minrowcor_main <- readRDS(paste0(simID, '--', sim, '/minrowcor_main.rds'))
    # BIC 1iter results
    BICscores_main1iter <- readRDS(paste0(simID, '--', sim, '/BICscores_1iter_main.rds'))  
    minrowcor_main1iter <- readRDS(paste0(simID, '--', sim, '/minrowcor_1iter_main.rds'))  
    # baseline 
    BICscores_baseline <- readRDS(paste0(simID, '--', sim, '/BICscores_baseline.rds')) 
    minrowcor_baseline <- readRDS(paste0(simID, '--', sim, '/minrowcor_baseline.rds'))
    # flipflop
    BICscores_ff <- readRDS(paste0(simID, '--', sim, '/BICscores_ff.rds'))
    minrowcor_ff <- readRDS(paste0(simID, '--', sim, '/minrowcor_ff.rds'))
    
    bestk_bic_main <- which.min(BICscores_main)
    bestk_cor_main <- which.min(minrowcor_main)
    bestk_bic_1iter <- which.min(BICscores_main1iter)
    bestk_cor_1iter <- which.min(minrowcor_main1iter)
    bestk_bic_baseline <- which.min(BICscores_baseline)
    bestk_cor_baseline <- which.min(minrowcor_baseline)
    bestk_bic_ff <- which.min(BICscores_ff)
    bestk_cor_ff <- which.min(minrowcor_ff)
    
    SHD_stats <- readRDS(paste0(simID, '--', sim, '/SHDclose.rds'))
    
    testX <- sim_X(
      vers = sim+100, 
      n = args$n,
      p = estimands$realp,
      omg.sq = estimands$omg.sq,
      sig = estimands$sig, 
      b = estimands$b
    )    
    
    testllres <- testll(
      testX$X, simID, sim,
      kmainbic = bestk_bic_main, 
      kmaincor = bestk_cor_main, 
      kbenchbic = bestk_bic_baseline, 
      kbenchcor = bestk_cor_baseline,
      k1iterbic = bestk_bic_1iter,
      k1itercor = bestk_cor_1iter,
      kffbic = bestk_bic_ff,
      kffcor = bestk_cor_ff
    )
    
    SHD_stats$shdXmain1iter['testll'] = testllres$bic_1iter
    SHD_stats$shdXmain1iterCor['testll'] = testllres$cor_1iter
    SHD_stats$shdXmain['testll'] = testllres$bic_main
    SHD_stats$shdXmainCor['testll'] = testllres$cor_main
    SHD_stats$shdXbaseline['testll'] = testllres$bic_baseline
    SHD_stats$shdXbaselineCor['testll'] = testllres$cor_baseline
    SHD_stats$shdXff['testll'] = testllres$bic_ff
    SHD_stats$shdXffCor['testll'] = testllres$cor_ff
    
    saveRDS(SHD_stats, file = paste0(simID, '--', sim, "/SHDclose.rds"))
  }
  setwd("~/Documents/research/dag_network")
}



process_output_unordered_tmp <- function(
  simID = simID,
  estimands,
  args,
  start=1,
  nsim=10, 
  thr = 0.1
){
  bstar_adj_cpdag <- bnstruct::dag.to.cpdag(1*(estimands$b != 0)) 
  p = dim(bstar_adj_cpdag)[1]
  dimnames(bstar_adj_cpdag) <- list(as.character(1:p), as.character(1:p))
  for(sim in start:nsim){
    setwd(paste0("output/",args$setting, "/", args$setting, "--", sim))
    cat(paste0('[INFO] Processing sim ', sim, '\n'))
    
    testX <- sim_X(
      vers = sim+100, 
      n = args$n,
      p = estimands$realp,
      omg.sq = estimands$omg.sq,
      sig = estimands$sig, 
      b = estimands$b
    ) 
    X = testX$X
    XXp <- Permute_X(X, seed = sim+1)
    Xp <- XXp$Xperm
    saveRDS(testX, file = "testX.rds")
    saveRDS(XXp, file = "testXp.rds")
    
    testllres <- testll_unordered(Xp)
    
    SHD_stats <- readRDS('SHDstats.rds')    
    SHD_stats$shd_pc['testll'] = testllres$pcll
    SHD_stats$shd_pc_decor['testll'] = testllres$pcdecorll
    SHD_stats$shd_ges['testll'] = testllres$gesll
    SHD_stats$shd_ges_decor['testll'] = testllres$gesdecorll
    SHD_stats$shd_sbn['testll'] = testllres$sbnll
    SHD_stats$shd_sbn_decor['testll'] = testllres$sbndecorll
    SHD_stats$shdXmain['testll'] = NA
    saveRDS(SHD_stats, file = 'SHDstats.rds')
    setwd("~/Documents/research/dag_network")  
  }
}


simIDs = c(
  '603','605'
)

# shd_average <- readRDS("~/Documents/research/dag_network/output/101/shd_average.rds")

setwd("~/Documents/research/dag_network")  
for(simID in simIDs){
  setwd(paste0("~/Documents/research/dag_network/output/", simID))
  args <- readRDS("args.rds")
  estimands <- readRDS("estimands.rds")
  # num_params = estimands$s0 + sum(estimands$theta!=0)
  testlldiffs = vector(mode = 'list', length = 10)  
  for(i in 1:10){
    SHDclose <- readRDS(paste0(simID, "--", i, "/SHDclose.rds"))
    testlldiffs[[i]] <- sapply(SHDclose, '[[', 'testll')
  }
  X <- readRDS(paste0(simID, "--", i, "/X.rds"))
  df <- data.frame(matrix(unlist(testlldiffs), nrow=length(testlldiffs), byrow=TRUE))
  colnames(df) <- names(SHDclose)
  dfplot = df %>% dplyr::select(shdXmain, shdXbaseline) %>% 
    gather(method, testll) %>% 
    mutate(testll = (testll - median(df$shdXbaseline))  / prod(dim(X$X))) 
  dfplot$method <- rep(c('BCD', 'Baseline'), each = 10)
  # dfplot$method <- rep(c('BCD', 'Baseline'), each = 10)
  plot_testll(dfplot, theta_type = simID)
  setwd("~/Documents/research/dag_network")  
}


for(simID in simIDs){
  shd_average <- readRDS(paste0("~/Documents/research/dag_network/output/", simID, "/shd_average.rds"))
  print(shd_average %>% dplyr::select(shdXmain, shdXmain1iter) %>% tail(2))
  print(simID)
}



# unordered cases ---------------------------------------------------------

n = 2
cols = gg_color_hue(n)

dev.new(width = 4, height = 4)
plot(1:n, pch = 16, cex = 2, col = cols)

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

simIDs = c(
  # '121001',
  # '121002',
  # '121003',
  # '121005',
# '121006')
# '701',
# '705',
'703'
)


setwd("~/Documents/research/dag_network")  
for(simID in simIDs){
  setwd(paste0("~/Documents/research/dag_network/output/", simID))
  args <- readRDS("args.rds")
  estimands <- readRDS("estimands.rds")
  # num_params = estimands$s0 + sum(estimands$theta!=0)
  testlldiffs = vector(mode = 'list', length = 10)  
  for(i in 1:10){
    SHDclose <- readRDS(paste0(simID, "--", i, "/SHDstats.rds"))
    testlldiffs[[i]] <- sapply(SHDclose, '[[', 'testll')
  }
  X <- readRDS(paste0(simID, "--", i, "/X.rds"))
  df <- data.frame(matrix(unlist(testlldiffs), nrow=length(testlldiffs), byrow=TRUE))
  colnames(df) <- names(SHDclose)
  dfplot = df %>% dplyr::select(-shdXmain) %>% 
    summarise(PC = shd_pc_decor - shd_pc,
              GES = shd_ges_decor - shd_ges,
              SBN = shd_sbn_decor - shd_sbn) %>% 
    gather(method, testll) %>% 
    mutate(testll = testll / prod(dim(X$X)))
  # dfplot$method <- rep(c('BCD', 'Baseline'), each = 10)
  # dfplot <- dfplot %>% filter(testll < 0.24, testll > 0.18)
  plot_testll(dfplot, theta_type = simID)
  setwd("~/Documents/research/dag_network")  
}

# theta error -------------------------------------------------------------
simIDs = c(
  '121001',
  '121002',
  '121003',
  '121005',
  '121006')
  # '101',
  # '002',
  # '003',
  # '007',
  # '004')


simIDs = c(
  '00011',
  '00022',
  '00033',
  '00055',
  '00044',
  '321001',
  '321002',
  '321003',
  '321004',
  '321005'
)
  


for(simID in simIDs){
  shd_average <- readRDS(paste0("~/Documents/research/dag_network/output/", simID, "/shd_average.rds"))
  print(shd_average %>% dplyr::select(shdXmain, shdXff, shdXmain1iter) %>% tail(2))
  print(simID)
}







