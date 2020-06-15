ThetaexpDecay <- genPartiallyConnect(n, nBlocks, bSizes, threashold = 0.05, seed = 12)
ThetaexpDecay$theta %>% cov2cor()
ThetaexpDecay$zeropos

test <- genTheta(10, graphType = 'genPartiallyConnect', threashold = 0.1, nBlocks = 2)
test

# # Network811 <- read.table("real_Sigma/Cross_Parker-Consulting_info.txt", quote = "\"", comment.char = "")
# freemans <- read.table("real_Sigma/Freemans_EIES-1_n48.txt", quote="\"", comment.char="")
# # USairport500 <- read.table("real_Sigma/USairport500.txt", quote="\"", comment.char="")


sample_corr <- runif(bSizes[4],-1,1)

mat <- matrix(c(0.5,-0.3,0.6,-0.3,0.5,0,0.6,0,0.5), 3, 3)
(mat + diag(3)*.4) %>% eigen() %>% "$"(values)

sig <- (mat + diag(3)*.4) %>% solve() 
sig

sig %>% solve() %>% round(10) 

sig %>% cov2cor() %>% solve() %>% round(10)

mat %>% solve()
eigen(mat)


cor1 <- diag(1 / sqrt(diag(sig)))
cor1
cor1inv <- diag(sqrt(diag(sig)))
cor1inv 


