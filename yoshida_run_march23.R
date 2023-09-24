#Packages
{
  library(Matrix) #for building random effects design matrix
  library(mvnfast) #for fast multivariate normal sampling
  library(MCMCpack) #for inverse-wishart sampling
  library(car) #for 95% prediction intervals
  library(plotrix) #for radial plots
  library(BayesLogit) #for Polya-Gamma auxiliary variables
}

#Internal algorithm functions
{
source("mmln_helpers.R")
}

#Main algorithms
{
source("mmln.R")
}

#The real data
#Updated every year, including all players to have played in both NPB and MLB
{
  
  dat <- read.csv("NZPAcrosslist_Final.csv")
  dat.test <- read.csv("NZPAcrosslist_Final_test.csv")
  dat <- dat[,-1]
  
}

#Prediction functions
{
source("mmln_predict.R")
}

#Data cleaning/preparation
{
  
  #The full, real data X matrix
  true_X <- dat[,c(3:5, 8:9, 11, 32, 33)]
  test_X <- dat.test[which(dat.test$switch == 0),c(4:6,9:10,12,33:34)]
  newp_X <- dat.test[which(dat.test$switch == 1),c(4:6,9:10,12,33:34)]
  
  #Batting Hand
  true_X$bats <- as.character(true_X$bats)
  test_X$bats <- as.character(test_X$bats)
  newp_X$bats <- as.character(newp_X$bats)
  true_X$bats[which(true_X$bats=="B")] <- as.character("S")
  test_X$bats[which(test_X$bats=="B")] <- as.character("S")
  newp_X$bats[which(newp_X$bats=="B")] <- as.character("S")
  true_X$bats[which(true_X$bats=="R")] <- as.character("H")
  test_X$bats[which(test_X$bats=="R")] <- as.character("H")
  newp_X$bats[which(newp_X$bats=="R")] <- as.character("H")
  true_X$bats <- as.numeric(as.factor(true_X$bats))
  test_X$bats <- as.numeric(as.factor(test_X$bats))
  newp_X$bats <- as.numeric(as.factor(newp_X$bats))
  
  #League
  true_X$lgID <- as.numeric(as.factor(true_X$lgID))
  test_X$lgID <- as.numeric(as.factor(test_X$lgID))
  newp_X$lgID <- as.numeric(as.factor(newp_X$lgID))
  
  #Position
  #Lump Infielders and Outfielders together, leave DH and P separate
  true_X$POS <- as.character(true_X$POS)
  true_X$POS[which(true_X$POS=="1B")] <- as.character("IF")
  true_X$POS[which(true_X$POS=="2B")] <- as.character("IF")
  true_X$POS[which(true_X$POS=="3B")] <- as.character("IF")
  true_X$POS[which(true_X$POS=="SS")] <- as.character("IF")
  true_X$POS[which(true_X$POS=="C")] <- as.character("IF")
  true_X$POS[which(true_X$POS=="CF")] <- as.character("OF")
  true_X$POS[which(true_X$POS=="LF")] <- as.character("OF")
  true_X$POS[which(true_X$POS=="RF")] <- as.character("OF")
  true_X$POS[which(true_X$POS=="P")] <- as.character("ZP")
  true_X$POS[which(true_X$POS=="DH")] <- as.character("XDH")
  true_X$POS <- as.numeric(as.factor(true_X$POS))
  
  test_X$POS <- as.character(test_X$POS)
  test_X$POS[which(test_X$POS=="1B")] <- as.character("IF")
  test_X$POS[which(test_X$POS=="2B")] <- as.character("IF")
  test_X$POS[which(test_X$POS=="3B")] <- as.character("IF")
  test_X$POS[which(test_X$POS=="SS")] <- as.character("IF")
  test_X$POS[which(test_X$POS=="C")] <- as.character("IF")
  test_X$POS[which(test_X$POS=="CF")] <- as.character("OF")
  test_X$POS[which(test_X$POS=="LF")] <- as.character("OF")
  test_X$POS[which(test_X$POS=="RF")] <- as.character("OF")
  test_X$POS[which(test_X$POS=="P")] <- as.character("ZP")
  test_X$POS[which(test_X$POS=="DH")] <- as.character("XDH")
  test_X$POS <- as.numeric(as.factor(test_X$POS))
  
  newp_X$POS <- as.character(newp_X$POS)
  newp_X$POS[which(newp_X$POS=="1B")] <- as.character("IF")
  newp_X$POS[which(newp_X$POS=="2B")] <- as.character("IF")
  newp_X$POS[which(newp_X$POS=="3B")] <- as.character("IF")
  newp_X$POS[which(newp_X$POS=="SS")] <- as.character("IF")
  newp_X$POS[which(newp_X$POS=="C")] <- as.character("IF")
  newp_X$POS[which(newp_X$POS=="CF")] <- as.character("OF")
  newp_X$POS[which(newp_X$POS=="LF")] <- as.character("OF")
  newp_X$POS[which(newp_X$POS=="RF")] <- as.character("OF")
  newp_X$POS[which(newp_X$POS=="P")] <- as.character("ZP")
  newp_X$POS[which(newp_X$POS=="DH")] <- as.character("XDH")
  newp_X$POS <- as.numeric(as.factor(newp_X$POS))
  
  #Create dummy columns
  A <- model.matrix(age ~ as.factor(bats), true_X)
  C <- model.matrix(age ~ as.factor(lgID), true_X)
  D <- model.matrix(age ~ as.factor(POS), true_X)
  true_X1 <- cbind(true_X[,c(1:2,4)], age2 = true_X[,4]^2, A[,2:3], C[,2:4], D[,2:4])
  colnames(true_X1)[5:12] <- c("bleft", "bswitch", "cl", "nl", "pl", "of", "dh", "p")
  
  At <- model.matrix(age ~ as.factor(bats), test_X)
  Ct <- model.matrix(age ~ as.factor(lgID), test_X)
  Dt <- model.matrix(age ~ as.factor(POS), test_X)
  test_X1 <- cbind(test_X[,c(1:2,4)], age2 = test_X[,4]^2, At[,2:3], Ct[,2:4], Dt[,2:4])
  rownames(test_X1) <- NULL
  colnames(test_X1)[5:12] <- c("bleft", "bswitch", "cl", "nl", "pl", "of", "dh", "p")
  
  An <- model.matrix(age ~ as.factor(bats), newp_X)
  Cn <- model.matrix(age ~ as.factor(lgID), newp_X)
  Dn <- model.matrix(age ~ as.factor(POS), newp_X)
  newp_X1 <- cbind(newp_X[,c(1:2,4)], age2 = newp_X[,4]^2, An[,2:3], Cn[,2:4], Dn[,2:4])
  rownames(newp_X1) <- NULL
  colnames(newp_X1)[5:12] <- c("bleft", "bswitch", "cl", "nl", "pl", "of", "dh", "p")
  
  #Remove added fake guys (will depend on the data set)
  test_X1 <- test_X1[-c(51:53),]
  newp_X1 <- newp_X1[-c(11:12),]
  
  #Scale
  scalecombo <- rbind(true_X1, test_X1, newp_X1)
  
  true_X2 <- cbind(scale(true_X1[,1:4]), true_X1[,5:12])
  #true_X2 <- cbind(scale(scalecombo[,1:4])[1:nrow(true_X1),], true_X1[,c(5:12)])
  head(true_X2)
  test_X2 <- cbind(scale(scalecombo[,1:4])[(nrow(true_X1)+1):(nrow(true_X1)+nrow(test_X1)),], test_X1[,c(5:12)])
  head(test_X2)
  newp_X2 <- cbind(scale(scalecombo[,1:4])[(nrow(true_X1)+nrow(test_X1)+1):nrow(scalecombo),], newp_X1[,c(5:12)])
  head(newp_X2)
  
  #Count matrix, W
  true_W <- as.matrix(cbind(dat[,c(17:20,25:27,29:30)], OIP = dat$PA - rowSums(dat[,c(17:20,25:27,29:30)])))
  head(true_W)
  test_W <- as.matrix(cbind(dat.test[which(dat.test$switch == 0),c(18:21,26:28,30:31)], OIP = dat.test$PA[which(dat.test$switch == 0)] - rowSums(dat.test[which(dat.test$switch == 0),c(18:21,26:28,30:31)])))
  rownames(test_W) <- NULL
  head(test_W)
  newp_W <- as.matrix(cbind(dat.test[which(dat.test$switch == 1),c(18:21,26:28,30:31)], OIP = dat.test$PA[which(dat.test$switch == 1)] - rowSums(dat.test[which(dat.test$switch == 1),c(18:21,26:28,30:31)])))
  rownames(newp_W) <- NULL
  head(newp_W)
  
  #Remove fake guys
  test_W <- test_W[-c(51:53),]
  newp_W <- newp_W[-c(11:12),]
  
  #Random effects design matrix, Z
  f <- factor(true_X$mergeon, levels = unique(true_X$mergeon))
  Ji <- t(as(f, Class = "sparseMatrix"))
  randX <- matrix(1, ncol = 1, nrow = nrow(true_X))
  Z <- t(KhatriRao(t(Ji), t(randX)))
  
  ft <- factor(test_X$mergeon, levels = unique(test_X$mergeon))
  Jit <- t(as(ft, Class = "sparseMatrix"))
  randXt <- matrix(1, ncol = 1, nrow = nrow(test_X))
  test_Z <- t(KhatriRao(t(Jit), t(randXt)))
  
  #Add the intercept
  true_X2 <- cbind(int = 1, true_X2, player = as.factor(as.numeric(factor(dat$mergeon, levels = unique(dat$mergeon)))))
  
  #mergeon without fakeos
  newmerge <- dat.test$mergeon[!c(dat.test$mergeon %in% c("guy 199", "guy2 199"))]
  newswitch <- dat.test$switch[!c(dat.test$mergeon %in% c("guy 199", "guy2 199"))]
  test_X2 <- cbind(int = 1, test_X2, player = as.factor(as.numeric(factor(newmerge[which(newswitch == 0)], levels = unique(newmerge))) + ncol(Z)))
  newp_X2 <- cbind(int = 1, newp_X2, player = as.factor(as.numeric(factor(newmerge[which(newswitch == 1)], levels = unique(newmerge))) + ncol(Z)))
  row.names(true_X2) <- NULL
  row.names(test_X2) <- NULL
  row.names(newp_X2) <- NULL
  
  
  pa <- rowSums(true_W)
  test_pa <- rowSums(test_W)
  newp_pa <- rowSums(newp_W)
  
}

#Run it (original prior)
#start.test <- Sys.time()
#testrun <- mixgibbs(W = true_W, X = as.matrix(true_X2[,-14]), Z = Z, base = 10, iter = 10000, proposal = "normbeta", whichYstart = "mean", whichCOVprior = "uninf")
#end.test <- Sys.time()

#Do it with Gelman's weakly informative prior
#Seems to converge much faster than original, completely disperse prior
#Still need to figure out if final prediction intervals are actually narrower though...
start.test <- Sys.time()
testrun <- mixgibbs(W = true_W, X = as.matrix(true_X2[,-14]), Z = Z, base = 10, iter = 100, proposal = "normbeta", whichYstart = "mean", whichCOVprior = "weak")
end.test <- Sys.time()

burnthin <- seq(301, 101, 2)

#Will use this run to predict the new players.
#To do so, need to get their random effects.

Y_samples <- testrun$Y[burnthin,]
Beta_samples <- testrun$Beta[burnthin,]
Sigma_samples <- testrun$Sigma[burnthin,]
Phi_samples <- testrun$Phi[burnthin,]
Psi_samples <- testrun$Psi[burnthin,]
MHRatio_samples <- testrun$MHRatio

#Psi (random effect) Recovery
#Careful, test_Z needs to be updated for fake guy(s)
test_psirecovery <- newobsrandeff(np_W = test_W, np_X = test_X2[,-14], np_Z = test_Z[1:nrow(test_W),-c(11:12)], randeff = 1, base = 10, iter = 100, mhiter = 1, mhscale = 1, proposal = "normbeta", sigma_chain = Sigma_samples, phi_chain = Phi_samples, beta_chain = Beta_samples, whichYstart = "mean", whichCOVprior = "weak", oldway = F)

#And posterior predictive draws
#Requires X matrix for new player seasons

#For in-sample players
truepredW <- true_W[which(dat$switch == 1),]
truepredX <- true_X2[which(dat$switch == 1),]
preds <- vector("list", nrow(truepredW))
for(i in 1:nrow(truepredW)){
  preds[[i]] <- predict.randeffseparate(np_pa = rowSums(truepredW)[i], np_X = t(truepredX[i,-14]), np_Z = matrix(1,1,1), dimW = 10, sigma_chain = Sigma_samples, phi_chain = Phi_samples, beta_chain = Beta_samples, psi_chain = Psi_samples[,c((9*i - 8):(9*i))], mean = TRUE)
}

#All seasons
fullpreds <- vector("list", nrow(true_W))
for(i in 1:nrow(true_W)){
  j <- as.numeric(true_X2$player)[i]
  fullpreds[[i]] <- predict.randeffseparate(np_pa = pa[i], np_X = t(true_X2[i,-14]), np_Z = matrix(1,1,1), dimW = 10, sigma_chain = Sigma_samples, phi_chain = Phi_samples, beta_chain = Beta_samples, psi_chain = Psi_samples[,c((9*j - 8):(9*j))], mean = TRUE)
}

# New Players
# Pre-season
true_newp_pa <- newp_pa
# After season
#true_newp_pa <- c(446, 484, 285, 465, 127) #update with new players after season ends

newp_W_news <- vector("list", nrow(newp_W))
for(i in 1:nrow(newp_W)){
  newp_W_news[[i]] <- predict.randeffseparate(np_pa = true_newp_pa[i], np_X = t(newp_X2[i,-14]), np_Z = matrix(1, 1, 1), dimW = 10, sigma_chain = Sigma_samples, phi_chain = Phi_samples, beta_chain = Beta_samples, psi_chain = test_psirecovery[,c((9*i - 8):(9*i))], mean = TRUE)
}

#How the predictions work
post.means <- matrix(unlist(lapply(preds, function(x) colMeans(x$predw))), byrow = T, ncol = 10)
colnames(post.means) <- colnames(truepredW)
head(post.means)

#Full
full.post.means <- matrix(unlist(lapply(fullpreds, function(x) colMeans(x$predw))), byrow = T, ncol = 10)
colnames(full.post.means) <- colnames(true_W)
head(full.post.means)

for(i in 1:10){
  plot(sqrt(true_W[,i]), sqrt(full.post.means[,i]))
}

for(i in 1:10){
  plot(sqrt(truepredW[,i]), sqrt(post.means[,i]))
}

for(i in 1:10){
  print(cor(sqrt(true_W[,i]), sqrt(full.post.means[,i])))
}

for(i in 1:10){
  print(cor(sqrt(truepredW[,i]), sqrt(post.means[,i])))
}

bavg <- function(w){
  H <- sum(w[1:4])
  AB <- sum(w[c(1:4,9:10)])
  return(H/AB)
}

obp <- function(w){
  AB <- sum(w[c(1:4,9:10)])
  OBPtop <- sum(w[c(1:4,7:8)])
  OBPbot <- AB + sum(w[c(6:8)])
  OBP <- OBPtop/OBPbot
  return(OBP)
}

slg <- function(w){
  AB <- sum(w[c(1:4,9:10)])
  TB <- sum(1:4*w[1:4])
  SLG <- TB/AB
  return(SLG)
}

ops <- function(w){
  AB <- sum(w[c(1:4,9:10)])
  OBPtop <- sum(w[c(1:4,7:8)])
  OBPbot <- AB + sum(w[c(6:8)])
  OBP <- OBPtop/OBPbot
  TB <- sum(1:4*w[1:4])
  SLG <- TB/AB
  OPS <- OBP + SLG
  return(OPS)
}

plot(apply(truepredW, 1, bavg), apply(post.means, 1, bavg))
plot(apply(truepredW, 1, obp), apply(post.means, 1, obp))
plot(apply(truepredW, 1, slg), apply(post.means, 1, slg))


ops.preds <- unlist(lapply(lapply(preds, function(x) apply(x$predw, 1, ops)), mean, na.rm = T))
plot(apply(truepredW[rowSums(truepredW) > 100,], 1, ops), ops.preds[rowSums(truepredW) > 100])
#plot(apply(truepredW, 1, ops), ops.preds)
abline(0,1)
cor(apply(truepredW[rowSums(truepredW) > 100,], 1, ops), ops.preds[rowSums(truepredW) > 100])
#cor(apply(truepredW, 1, ops), ops.preds)

plot(apply(true_W[pa > 500,], 1, ops), apply(full.post.means[pa > 500,], 1, ops))
abline(0,1)
cor(apply(true_W[pa > 500,], 1, ops), apply(full.post.means[pa > 500,], 1, ops))


#Katoh
round(colMeans(newp_W_news[[1]]$predw), 1)
bavg(colMeans(newp_W_news[[1]]$predw))
obp(colMeans(newp_W_news[[1]]$predw))
slg(colMeans(newp_W_news[[1]]$predw))
ops(colMeans(newp_W_news[[1]]$predw))
#true_katoh <- c(66, 22, 2, 14, 0, 3, 42, 4, 110, 183) #update with new data after season ends
#bavg(true_katoh)

#Yoshida
round(colMeans(newp_W_news[[2]]$predw), 1)
bavg(colMeans(newp_W_news[[2]]$predw))
obp(colMeans(newp_W_news[[2]]$predw))
slg(colMeans(newp_W_news[[2]]$predw))
ops(colMeans(newp_W_news[[2]]$predw))
#true_yoshida <- c(59, 21, 1, 24, 0, 3, 42, 1, 109, 224) #update with new data after season ends

#Schwindel
round(colMeans(newp_W_news[[3]]$predw), 1)
bavg(colMeans(newp_W_news[[3]]$predw))
obp(colMeans(newp_W_news[[3]]$predw))
slg(colMeans(newp_W_news[[3]]$predw))
ops(colMeans(newp_W_news[[3]]$predw))
#true_alcantara <- c(35, 6, 0, 14, 0, 1, 21, 0, 87, 120) #update with new data after season ends

#Brinson
round(colMeans(newp_W_news[[4]]$predw), 1)
bavg(colMeans(newp_W_news[[4]]$predw))
obp(colMeans(newp_W_news[[4]]$predw))
slg(colMeans(newp_W_news[[4]]$predw))
ops(colMeans(newp_W_news[[4]]$predw))
#true_ogrady <- c(47, 24, 0, 15, 0, 2, 54, 5, 129, 189) #update with new data after season ends

#Gonzalez
round(colMeans(newp_W_news[[5]]$predw), 1)
bavg(colMeans(newp_W_news[[5]]$predw))
obp(colMeans(newp_W_news[[5]]$predw))
slg(colMeans(newp_W_news[[5]]$predw))
ops(colMeans(newp_W_news[[5]]$predw))
#true_valera <- c(17, 5, 0, 1, 1, 0, 13, 1, 11, 78) #update with new data after season ends

#Payton
round(colMeans(newp_W_news[[6]]$predw), 1)
bavg(colMeans(newp_W_news[[6]]$predw))
obp(colMeans(newp_W_news[[6]]$predw))
slg(colMeans(newp_W_news[[6]]$predw))
ops(colMeans(newp_W_news[[6]]$predw))
#true_valera <- c(17, 5, 0, 1, 1, 0, 13, 1, 11, 78) #update with new data after season ends

#Aquino
round(colMeans(newp_W_news[[7]]$predw), 1)
bavg(colMeans(newp_W_news[[7]]$predw))
obp(colMeans(newp_W_news[[7]]$predw))
slg(colMeans(newp_W_news[[7]]$predw))
ops(colMeans(newp_W_news[[7]]$predw))
#true_valera <- c(17, 5, 0, 1, 1, 0, 13, 1, 11, 78) #update with new data after season ends

#Davidson
round(colMeans(newp_W_news[[8]]$predw), 1)
bavg(colMeans(newp_W_news[[8]]$predw))
obp(colMeans(newp_W_news[[8]]$predw))
slg(colMeans(newp_W_news[[8]]$predw))
ops(colMeans(newp_W_news[[8]]$predw))
#true_valera <- c(17, 5, 0, 1, 1, 0, 13, 1, 11, 78) #update with new data after season ends

#Astudillo
round(colMeans(newp_W_news[[9]]$predw), 1)
bavg(colMeans(newp_W_news[[9]]$predw))
obp(colMeans(newp_W_news[[9]]$predw))
slg(colMeans(newp_W_news[[9]]$predw))
ops(colMeans(newp_W_news[[9]]$predw))
#true_valera <- c(17, 5, 0, 1, 1, 0, 13, 1, 11, 78) #update with new data after season ends

#MacKinnon
round(colMeans(newp_W_news[[10]]$predw), 1)
bavg(colMeans(newp_W_news[[10]]$predw))
obp(colMeans(newp_W_news[[10]]$predw))
slg(colMeans(newp_W_news[[10]]$predw))
ops(colMeans(newp_W_news[[10]]$predw))
#true_valera <- c(17, 5, 0, 1, 1, 0, 13, 1, 11, 78) #update with new data after season ends



#Credible Interval for Mean Prediction
ppsid_like <- function(ppsid, mu, phi){
  like <- dmvn(ppsid, mu = mu, sigma = phi, log = TRUE)
  return(like)
}

PostPhi <- matrix(colMeans(Phi_samples), byrow = T, ncol = 9)
PostBeta <- matrix(colMeans(Beta_samples), byrow = T, ncol = 9)
#Do it for Yoshida
meanpostpsi <- colMeans(newp_W_news[[2]]$psisamples)
LIKES <- apply(newp_W_news[[2]]$psisamples, 1, ppsid_like, mu = meanpostpsi, PostPhi)
pct95 <- round(.95*nrow(newp_W_news[[2]]$psisamples))
sortLIKES <- newp_W_news[[2]]$psisamples[order(LIKES, decreasing = T),]
predpsi95like <- sortLIKES[1:pct95,]
predXBpP <- matrix(0, ncol = 9, nrow = pct95)
for(i in 1:pct95){
  predXBpP[i,] <- matrix(unlist(newp_X2[2,-14]), nrow = 1)%*%PostBeta + predpsi95like[i,]
}
predpsi95prob <- cbind(exp(predXBpP)/(1 + rowSums(exp(predXBpP))), 1/(1 + rowSums(exp(predXBpP))))
countsint <- 500*predpsi95prob #at end of season, update 500 with true PA
round(apply(countsint, 2, range))
range(apply(countsint, 1, bavg))
range(apply(countsint, 1, obp))
range(apply(countsint, 1, slg))
range(apply(countsint, 1, ops))

# Season isn't over yet; once it is, update true_yoshida above and create the plots
df_newp <- tibble(
  stat = c(colnames(test_W), 'AVG', 'OBP', 'SLG', 'OPS'),
  truth =  c(true_yoshida, bavg(true_yoshida), obp(true_yoshida), slg(true_yoshida), ops(true_yoshida)),
  pred = c(round(colMeans(newp_W_news[[2]]$predw), 1), bavg(colMeans(newp_W_news[[2]]$predw)), obp(colMeans(newp_W_news[[2]]$predw)), slg(colMeans(newp_W_news[[2]]$predw)), ops(colMeans(newp_W_news[[2]]$predw))),
  pred_lower = c(round(apply(countsint, 2, range))[1,], range(apply(countsint, 1, bavg))[1], range(apply(countsint, 1, obp))[1], range(apply(countsint, 1, slg))[1], range(apply(countsint, 1, ops))[1]),
  pred_upper = c(round(apply(countsint, 2, range))[2,], range(apply(countsint, 1, bavg))[2], range(apply(countsint, 1, obp))[2], range(apply(countsint, 1, slg))[2], range(apply(countsint, 1, ops))[2])
)

df_newp[1:10,] %>%
  ggplot(aes(stat, pred)) +
  geom_point(col = 'deepskyblue4') +
  geom_errorbar(aes(ymin = pred_lower, ymax = pred_upper), width = 0.3) +
  geom_point(aes(y=truth),colour="red",pch='x',size = 7, alpha = .7) + 
  ggtitle('Masataka Yoshida Predictions vs. Truth 2023', subtitle = 'Using Mixed-MLN Model (Gerber and Craig, 2021)') +
  xlab('Statistic') +
  ylab('Value') +
  theme_grey(base_size = 16)

df_newp[11:14,] %>%
  ggplot(aes(stat, pred)) +
  geom_point(col = 'deepskyblue4', size = 2) +
  geom_errorbar(aes(ymin = pred_lower, ymax = pred_upper), width = 0.3) +
  geom_point(aes(y=truth),colour="red",pch='x',size = 7, alpha = .7) + 
  ggtitle('Masataka Yoshida Predictions vs. Truth 2023', subtitle = 'Using Mixed-MLN Model (Gerber and Craig, 2021)') +
  xlab('Statistic') +
  ylab('Value') +
  theme_grey(base_size = 16)
