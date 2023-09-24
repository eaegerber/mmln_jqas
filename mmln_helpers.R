#Helper functions for mmln

#Starting Values of Y Function
initialY <- function(W, base, zedsc = .05){
  Wmat <- as.matrix(W)
  N <- nrow(Wmat) #number of observations
  K <- ncol(Wmat) #number of multinomial categories
  M <- apply(Wmat, 1, sum) #exposure (plate appearances)
  
  tempZmat <- matrix(0,N,K)
  Ymat <- matrix(0,N,K-1)
  
  #dealing with zero counts
  for (i in 1:N){
    zeros <- which(Wmat[i,]==0)
    nzeros <- which(Wmat[i,]!=0)
    tempZmat[i,zeros] <- (Wmat[i,zeros]+zedsc)/(M[i]+zedsc*length(zeros))
    tempZmat[i,nzeros] <- (Wmat[i,nzeros])/(M[i]+zedsc*length(zeros))
    Ymat[i,] <- log(tempZmat[i,-base]/tempZmat[i,base])
  }
  
  return(Ymat)
}

newinitialY <- function(W, base){
  Wmat <- as.matrix(W)
  N <- nrow(Wmat)
  K <- ncol(Wmat)
  M <- apply(Wmat, 1, sum)
  
  tempZmat <- Wmat/M
  Ymat <- matrix(0, nrow = N, ncol = K-1)
  
  for(i in 1:N){
    if(length(which(tempZmat[i,]==0)) > 0){
      tempZmat[i,] <- (tempZmat[i,]*(M[i] - 1) + (1/K)) / M[i]
    }
    Ymat[i,] <- log(tempZmat[i,-base]/tempZmat[i,base])
  }
  
  return(Ymat)
  
}

#Inverse of matrix using Cholesky Decomposition (faster than solve())
cholinv <- function(x){
  chol2inv(chol(x))
}

#For the Static Beta Proposal
pstartoy <- function(pstarvec){
  y <- numeric(length=length(pstarvec))
  for(i in 1:length(pstarvec)){
    y[i] <- log(pstarvec[i]/(1-pstarvec[i]))
  }
  return(y)
}

ytopstar <- function(yvec){
  pstar <- numeric(length=length(yvec))
  for(i in 1:length(yvec)){
    pstar[i] <- exp(yvec[i])/(1 + exp(yvec[i]))
  }
  return(pstar)
}

betapropdist <- function(WMu_vec, Sigma){
  k <- ncol(Sigma)
  w_vec <- WMu_vec[1:(k+1)]
  mu_vec <- WMu_vec[(k+2):length(WMu_vec)]
  result <- numeric(length=k)
  for(i in 1:k){
    alpha <- ((1 + exp(mu_vec[i])) / Sigma[i,i]) - (exp(mu_vec[i]) / (1 + exp(mu_vec[i])))
    if(alpha < 0){
      alpha <- alpha + 1
    }
    beta <- alpha*exp(-mu_vec[i])
    alpha_star <- w_vec[i] + alpha
    beta_star <- w_vec[k+1] + beta
    result[i] <- rbeta(1, alpha_star, beta_star)
  }
  return(result)
}

betaloglike <- function(WMuPstar_vec, Sigma){
  k <- ncol(Sigma)
  w_vec <- WMuPstar_vec[1:(k+1)]
  mu_vec <- WMuPstar_vec[(k+2):(2*k+1)]
  pstar_vec <- WMuPstar_vec[(2*k+2):length(WMuPstar_vec)]
  loglike <- 0
  logjac <- 0
  for(i in 1:k){
    alpha <- ((1 + exp(mu_vec[i])) / Sigma[i,i]) - (exp(mu_vec[i]) / (1 + exp(mu_vec[i])))
    if(alpha < 0){
      alpha <- alpha + 1
    }
    beta <- alpha*exp(-mu_vec[i])
    alpha_star <- w_vec[i] + alpha
    beta_star <- w_vec[k+1] + beta
    loglike <- loglike + dbeta(pstar_vec[i], alpha_star, beta_star, log = T)
    logjac <- logjac + log(pstar_vec[i]) + log(1 - pstar_vec[i])
  }
  return(loglike + logjac)
}

#For the Static Normal Approx to Beta Proposal
normbetapropdist <- function(WMu_vec, Sigma){
  k <- ncol(Sigma)
  w_vec <- WMu_vec[1:(k+1)]
  mu_vec <- WMu_vec[(k+2):length(WMu_vec)]
  result <- numeric(length=k)
  for(i in 1:k){
    alpha <- ((1 + exp(mu_vec[i])) / Sigma[i,i]) - (exp(mu_vec[i]) / (1 + exp(mu_vec[i])))
    beta <- alpha*exp(-mu_vec[i])
    alpha_star <- w_vec[i] + alpha
    beta_star <- w_vec[k+1] + beta
    muprop <- digamma(alpha_star) - digamma(beta_star)
    sigprop <- sqrt(trigamma(alpha_star) + trigamma(beta_star))
    result[i] <- rnorm(1, muprop, sigprop)
  }
  return(result)
}

normbetaloglike <- function(WMuY_vec, Sigma){
  k <- ncol(Sigma)
  w_vec <- WMuY_vec[1:(k+1)]
  mu_vec <- WMuY_vec[(k+2):(2*k+1)]
  y_vec <- WMuY_vec[(2*k+2):length(WMuY_vec)]
  result <- 0
  for(i in 1:k){
    alpha <- ((1 + exp(mu_vec[i])) / Sigma[i,i]) - (exp(mu_vec[i]) / (1 + exp(mu_vec[i])))
    beta <- alpha*exp(-mu_vec[i])
    alpha_star <- w_vec[i] + alpha
    beta_star <- w_vec[k+1] + beta
    muprop <- digamma(alpha_star) - digamma(beta_star)
    sigprop <- sqrt(trigamma(alpha_star) + trigamma(beta_star))
    result <- result + dnorm(y_vec[i], muprop, sigprop, log = T)
  }
  return(result)
}

#Quickly convert Y to Pi
ytopi <- function(vec){
  c(exp(vec)/(1 + sum(exp(vec))), 1/(1 + sum(exp(vec))))
}