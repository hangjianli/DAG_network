# flipflop ----------------------------------------------------------------
# set.seed(1)
# resultflip <- flipflop(X = XX, n = args$n, p = estimands$realp,
#                        zeropos_theta = estimands$zeropos,
#                        max.iter = 20,
#                        lambda = lambda,
#                        lambda2 = 0.01)
# #
# psi_hat <- resultflip$psi_est
# theta_hat <- resultflip$theta_est
# psi_zeros <- which(abs(psi_hat) < 1e-4, arr.ind = T)
# theta_zeros <- which(abs(theta_hat) < 1e-4, arr.ind = T)
# if(any(diff(resultflip$likeli_seq[resultflip$likeli_seq!=0]) > 0))
#   warning("likelihood not descreasing ! \n")
# 
# png(paste0("B_est_flipflop_lik_", k ,".png"))
# plot(resultflip$likeli_seq[resultflip$likeli_seq!=0], pch = 16, type = "b")
# lines(resultflip$likeli_seq[resultflip$likeli_seq!=0], type = "o")
# dev.off()
# #
# cat("[INFO] fliplop MLE calculation. \n")
# resultMLE <- flipflop(X=XX, n = n, p = p,
#                       zeropos_theta = theta_zeros,
#                       zeropos_psi = psi_zeros,
#                       max.iter = 15,
#                       lambda = 0.0001,
#                       lambda2 = 0.0001)
# psi_mle <- resultMLE$psi_est
# b_est <- resultMLE$b_est
# theta_mle <- resultMLE$theta_est
# # 
# heatmap.2(abs(solve(theta_mle)),
#           dendrogram = "none",
#           Rowv = F, Colv = F,
#           main= paste0("flipflop_", k),
#           trace = "none")
# # 
# saveRDS(resultMLE, file = paste0(format(today_date, "%Y-%m-%d"),"-lam-",
#                                  k,"-flipflopMLE",".rds"))
# BICscores_flipflop[k] <- -n*log(det(psi_mle)) - p*log(det(theta_mle)) +
#   sum(diag(theta_mle%*%XX%*%psi_mle%*%t(XX))) + log(max(n,p))*sum(abs(psi_mle) > 1e-4)
# cor_est <- cor(t(chol(resultMLE$theta_est)%*%(XX - XX%*%b_est)))
# # sd_fliplfop[k] <- sd(cor_est[upper.tri(cor_est)])
# minrowcor_flip[k] <- sum(abs(cor_est[upper.tri(cor_est)]))
# saveRDS(minrowcor_flip,"minrowcor_flip.rds")
# saveRDS(BICscores_flipflop,"BICscores_flipflop.rds")