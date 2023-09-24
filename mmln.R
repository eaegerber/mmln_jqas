#Main mmln functions including likelihood tracking

#W is the N by (d + 1) Multinomial Count Matrix
#X is the N by p Covariate Matrix with Intercept
#Z is the N by m*q Random Effects Design Matrix
#randeff = q, the number of random effects
#base is which of the (d + 1) categories to use as a baseline
#iter is the total number of Gibbs iterations
#mhiter is the number of Metropolis-Hastings iterations for sampling Y
#mhscale scales the proposal variance in the MH in the random walk case
#proposal determines which proposal scheme to use
#block determines whether to use the block update or individual gibbs steps for location parameters
#r.seed sets the seed for the Gibbs run
#trY,...,trPhi are the true values of the respective parameters (troubleshooting in simulation studies)
#startru determines whether to start parameters at their true values (troubleshooting; requires true values)
#fixY,...,fixPhi are logicals fixing the respective parameters at their true values for the duration of the algorithm (troubleshooting; requires true values)
#test pulls you into the function at a chosen place (troubleshooting)
#setting mixed = FALSE will make it a Fixed effects model (no random effects). It will still require a Z matrix (if I default it to NULL than the default is not mixed).

mixgibbs <- function(W, X, Z, randeff = 1, base, iter = 1000, mhiter = 1, mhscale = 1, proposal = c("norm", "beta", "normbeta"), block = TRUE, r.seed = 42, trY = NULL, trBeta = NULL, trPsi = NULL, trSigma = NULL, trPhi = NULL, startru = FALSE, fixY = FALSE, fixBeta = FALSE, fixPsi = FALSE, fixSigma = FALSE, fixPhi = FALSE, whichYstart = c("mean", "inity"), whichCOVprior = c("uninf", "weak"), mixed = TRUE, test = FALSE){
  
  #Set seed for replicability
  set.seed(r.seed)
  
  #Mixed Model
  if(mixed == TRUE){
    
    #This if statement can be moved around for troubleshooting at different points
    if(test == TRUE){
      #if(i == 1){
      browser()         
      #}
    }
    
    #Make sure base is the last category
    W <- cbind(W[,-base], W[,base])
    base <- ncol(W)
    
    #Some dimensions
    p <- ncol(X) #number of covariates, including intercept
    d <- ncol(W) - 1 #number of multivariate categories
    N <- nrow(X) #number of total observations
    q <- randeff #number of random player effects (default is simply an intercept)
    m <- ncol(Z)/randeff #number of players
    
    #Algorithm uses a wide X
    X <- t(X) #should be p by N
    
    #Get initial estimates of the log-odds
    if(fixY == FALSE){
      Y <- t(newinitialY(W, base)) #should be d by N
    } else{
      Y <- t(trY)
    }
    
    #Beta set-up
    if(fixBeta == FALSE){
      prior_beta <- t(tcrossprod(cholinv(tcrossprod(X)), tcrossprod(Y, X)))
      betatrack <- matrix(0, nrow = iter + 1, ncol = d*p)
      betatrack[1,] <- as.vector(prior_beta)
      if(startru == TRUE){
        prior_beta <- trBeta
        betatrack <- matrix(0, nrow = iter + 1, ncol = d*p)
        betatrack[1,] <- as.vector(prior_beta)
      }
    } else{
      prior_beta <- trBeta
      betatrack <- matrix(0, nrow = iter + 1, ncol = d*p)
      betatrack[1,] <- as.vector(prior_beta)
    }
    
    #Psi set-up
    if(fixPsi == FALSE){
      #prior_psi <- apply(rbind(colSums(Z), c(0, colSums(Z)[-m])), 2, function(x) colMeans(t(Y - tcrossprod(prior_beta, t(X)))[(x[2]+1):(x[1]+x[2]),]) )
      prior_psi <- matrix(0, nrow = d, ncol = q*m)
      psitrack <- matrix(0, nrow = iter + 1, ncol = length(prior_psi))
      psitrack[1,] <- as.vector(prior_psi)
      if(startru == TRUE){
        prior_psi <- t(trPsi) #should be d by q*m
        psitrack <- matrix(0, nrow = iter + 1, ncol = length(prior_psi))
        psitrack[1,] <- as.vector(prior_psi)
      }
    } else{
      prior_psi <- t(trPsi)
      psitrack <- matrix(0, nrow = iter + 1, ncol = length(prior_psi))
      psitrack[1,] <- as.vector(prior_psi)
    }
    
    #Sigma set-up
    #Hyperparameters for prior distribution of Sigma
    if(whichCOVprior == "weak"){
      v_1 <- d + 2 #Gelman's weakly informative prior
      Lambda_1 <- diag(d)*(1/(2*.001)) #Gelman's weakly informative prior
    }else{
      v_1 <- d #completely uninformatve prior
      Lambda_1 <- diag(d) + matrix(1, nrow = d)%*%matrix(1, ncol = d) #completely uninformatve prior
      
    }
    
    if(fixSigma == FALSE){
      prior_sigma <- diag(d)
      sigmatrack <- matrix(0, nrow = iter + 1, ncol = length(prior_sigma))
      sigmatrack[1,] <- as.vector(prior_sigma)
      if(startru == TRUE){
        prior_sigma <- trSigma
        sigmatrack <- matrix(0, nrow = iter + 1, ncol = length(prior_sigma))
        sigmatrack[1,] <- as.vector(prior_sigma)
      }
    } else{
      prior_sigma <- trSigma
      sigmatrack <- matrix(0, nrow = iter + 1, ncol = length(prior_sigma))
      sigmatrack[1,] <- as.vector(prior_sigma)
    }
    
    #Phi set-up
    #Hyperparameters for prior distribution of Phi
    if(whichCOVprior == "weak"){
      v_2 <- d*q + 2 #Gelman's weakly informative prior
      Lambda_2 <- diag(d*q)*(1/(2*.001)) #Gelman's weakly informative prior
    }else{
      v_2 <- d*q #completely uninformatve prior
      Lambda_2 <- diag(d*q) + matrix(1, nrow = d*q)%*%matrix(1, ncol = d*q)  #completely uninformatve prior
    }
    
    if(fixPhi == FALSE){
      prior_phi <- diag(d*q)
      phitrack <- matrix(0, nrow = iter + 1, ncol = length(prior_phi))
      phitrack[1,] <- as.vector(prior_phi)
      if(startru == TRUE){
        prior_phi <- trPhi
        phitrack <- matrix(0, nrow = iter + 1, ncol = length(prior_phi))
        phitrack[1,] <- as.vector(prior_phi)
      }
    } else{
      prior_phi <- trPhi
      phitrack <- matrix(0, nrow = iter + 1, ncol = length(prior_phi))
      phitrack[1,] <- as.vector(prior_phi)
    }
    
    #Get initial log-odds means, and then set starting values of Y to that (except doing that leads to horrible MH acceptance ratio when moving to more categories)
    prior_mu <- tcrossprod(prior_beta, t(X)) + prior_psi%*%as.matrix(t(Z))
    
    #Starting values for Y
    if(fixY == FALSE){
      if(whichYstart == "mean"){
        Y <- prior_mu
      }
    } else{
      Y <- t(trY)
    }
    
    ytrack <- matrix(0, nrow = iter + 1, ncol = d*N)
    ytrack[1,] <- as.vector(Y)
    
    #Certain values will stay fixed throughout the Gibbs chain
    v_1n <- v_1 + N - d*p #posterior df of Sigma
    v_2n <- v_2 + m #posterior df of Phi
    S_xx <- tcrossprod(X) #X sums of squares matrix
    
    #Finally, also want to keep track of MH acceptance probabilities for proposing Y
    MHacceptancelist <- matrix(0, nrow = iter, ncol = mhiter)
    
    #Full posterior log-likelihood and error term tracking can be recovered from the resulting output matrices, and to save computation time are ommitted from the main function
    
    #The MH-within-Gibbs Algorithm
    for(i in 1:iter){
      
      #        if(i %in% seq(100, iter, round(iter/100))){
      #          print(i)
      #        }
      
      #Metropolis-Hastings step
      #Y|W, X, Z, Beta, Psi, Sigma, Phi
      
      if(fixY == FALSE){
        for(j in 1:mhiter){
          Y <- t(Y) #put it long ways for proposing
          Mu <- t(prior_mu) #put it long ways for proposing
          wmu <- cbind(W, Mu) #for the proposal function
          wmuy <- cbind(wmu, Y)
          
          if(proposal == "norm"){
            Y_new <- Y + rmvn(N, integer(length = d), mhscale*prior_sigma) #propose centered at old Y
          }
          
          if(proposal == "beta"){
            Pstar <- t(apply(X = Y, MARGIN = 1, FUN = ytopstar)) #old Pstar
            Pstar_new <- t(apply(X = wmu, MARGIN = 1, FUN = betapropdist, Sigma = mhscale*prior_sigma))
            Y_new <- t(apply(X = Pstar_new, MARGIN = 1, FUN = pstartoy)) #newY based on proposed probabilities
            wmuy_new <- cbind(wmu, Y_new)
            wmuPstar <- cbind(wmu, Pstar)
            wmuPstar_new <- cbind(wmu, Pstar_new)
          }
          
          if(proposal == "normbeta"){
            Y_new <- t(apply(X = wmu, MARGIN = 1, FUN = normbetapropdist, Sigma = mhscale*prior_sigma)) #Fixed Normal proposal based on logit(Beta) expected values and variances
            wmuy_new <- cbind(wmu, Y_new)
          }
          
          sexpY <- rowSums(exp(Y))
          sexpYn <- rowSums(exp(Y_new))
          
          newnormpart <- .5*tcrossprod(tcrossprod((Y_new - Mu) , t(cholinv(prior_sigma))), (Y_new - Mu))
          oldnormpart <- .5*tcrossprod(tcrossprod((Y - Mu) , t(cholinv(prior_sigma))), (Y - Mu))
          newloglike <- rowSums(as.matrix(W[,1:d]*Y_new[,1:d])) - rowSums(W*(log1p(sexpYn))) - newnormpart[cbind(1:nrow(newnormpart), 1:nrow(newnormpart))]
          oldloglike <- rowSums(as.matrix(W[,1:d]*Y[,1:d])) - rowSums(W*(log1p(sexpY))) - oldnormpart[cbind(1:nrow(oldnormpart), 1:nrow(oldnormpart))]
          
          rm(newnormpart, oldnormpart) #remove from memory to help efficiency
          
          if(proposal == "norm"){
            #Random walk is symmetric, Metropolis, not MH
            logptolike <- 0
            logotplike <- 0
            
            rm(wmu, wmuy)
          }
          
          if(proposal == "beta"){
            logptolike <- apply(X = wmuPstar, MARGIN = 1, FUN = betaloglike, Sigma = mhscale*prior_sigma) #going to old p_star under proposal distribution, including Jacobian
            logotplike <- apply(X = wmuPstar_new, MARGIN = 1, FUN = betaloglike, Sigma = mhscale*prior_sigma) #going to proposed p_star under proposal, including Jacobian
            
            rm(wmu, wmuy, wmuy_new, wmuPstar, wmuPstar_new)
          }
          
          if(proposal == "normbeta"){
            logptolike <- apply(X = wmuy, MARGIN = 1, FUN = normbetaloglike, Sigma = mhscale*prior_sigma)
            logotplike <- apply(X = wmuy_new, MARGIN = 1, FUN = normbetaloglike, Sigma = mhscale*prior_sigma)
            
            rm(wmu, wmuy, wmuy_new)
          }
          
          ratio <- newloglike - oldloglike + logptolike - logotplike
          temp <- runif(N)
          difr <- as.numeric((ratio - log(temp)) >= 0)
          MHacceptancelist[i,j] <- mean(difr)
          Y <- difr*Y_new + (1 - difr)*Y #accept or reject samples to get draw
          Y <- t(Y) #put it back for the rest of the Gibbs
        }
        
        ytrack[i+1,] <- as.vector(Y)
      } else{
        #To Troubleshoot: keep Y fixed
        ytrack[i+1,] <- as.vector(Y)
      }
      
      #Using sampled Y, create necessary sums of squares matrices and OLS mean estimates for:
      #individual update of Beta
      if(fixBeta == FALSE){
        S_yx <- tcrossprod((Y - prior_psi%*%as.matrix(t(Z))), X)
        beta_hat <- tcrossprod(S_yx, t(cholinv(S_xx)))
        
        #block update of Beta
        margS_yx <- tcrossprod(Y, X)
        margbeta_hat <- tcrossprod(margS_yx, t(cholinv(S_xx)))
      }
      
      #Block update of Beta and Psi is more efficient than individual updates
      
      if(block == TRUE){
        
        #Marginal Update of Beta, not conditional on Psi
        #Beta|Y, X, Z, Sigma, Phi
        if(fixBeta == FALSE){
          sig_beta_inv <- matrix(0, ncol = p*d, nrow = p*d)
          for(k1 in 1:m){
            Zi <- as.matrix(t(Z[which(t(Z)[(k1*q - (q - 1)),]==1),(k1*q - (q - 1)):(k1*q)]))
            n_i <- ncol(Zi)
            Xi <- X[,which(t(Z)[(k1*q - (q - 1)),]==1)]
            if(nrow(X) == 1){
              Xi <- matrix(Xi, nrow = 1)
            }
            sig_beta_inv <- sig_beta_inv + kronecker(tcrossprod(Xi), cholinv(n_i*prior_phi + prior_sigma))
          }
          
          sig_beta <- cholinv(sig_beta_inv)
          mu_beta <- margbeta_hat
          
          post_beta_vec <- rmvn(1, mu_beta, sig_beta)
          prior_beta <- matrix(post_beta_vec, nrow = d) #should be d by p
          betatrack[i+1,] <- as.vector(prior_beta)
        } else{
          #To Troubleshoot: keep Beta fixed
          post_beta <- prior_beta
          betatrack[i+1,] <- as.vector(prior_beta)
        }
        
        #Psi| Beta, Y, X, Z, Sigma, Phi
        if(fixPsi == FALSE){
          post_psi <- matrix(0, nrow = d*q, ncol = m)
          for(k in 1:m){
            Xi <- X[,which(t(Z)[(k*q - (q - 1)),]==1)]
            Yi <- Y[,which(t(Z)[(k*q - (q - 1)),]==1)]
            Zi <- as.matrix(t(Z[which(t(Z)[(k*q - (q - 1)),]==1),(k*q - (q - 1)):(k*q)]))
            Ui <- cholinv(cholinv(prior_phi) + kronecker(cholinv(prior_sigma) , as.matrix(Zi%*%t(Zi))))
            
            sigZi <- kronecker(cholinv(prior_sigma), Zi)
            errvec <- matrix(as.vector(t(Yi - prior_beta%*%Xi)), byrow = T, ncol = 1)
            psi_hat <- Ui%*%sigZi%*%errvec
            
            post_psivec <- rmvn(1, as.vector(psi_hat), Ui)
            post_psi[,k] <- post_psivec
          }
          prior_psi_star <- post_psi #need it in alternate form for sampling Phi (d*q by m)
          post_psi <- matrix(as.vector(post_psi), nrow = d) #for the rest of the time (d by q*m)
          psitrack[i+1,] <- as.vector(post_psi)
          prior_psi <- post_psi
        } else{
          #To troubleshoot: keep Psi fixed
          prior_psi_star <- matrix(as.vector(prior_psi), ncol = m) #Not entirely sure this will work with multiple random effects
          post_psi <- prior_psi
          psitrack[i+1,] <- as.vector(post_psi)
        }
        
      }
      
      #Individual updates of Beta and Psi change the order a bit, to improve efficiency
      #Psi| Beta, Y, X, Z, Sigma, Phi
      
      if(block == FALSE){
        
        #Same as in block update
        if(fixPsi == FALSE){
          post_psi <- matrix(0, nrow = d*q, ncol = m)
          for(k in 1:m){
            Xi <- X[,which(t(Z)[(k*q - (q - 1)),]==1)]
            Yi <- Y[,which(t(Z)[(k*q - (q - 1)),]==1)]
            Zi <- as.matrix(t(Z[which(t(Z)[(k*q - (q - 1)),]==1),(k*q - (q - 1)):(k*q)]))
            Ui <- cholinv(cholinv(prior_phi) + kronecker(cholinv(prior_sigma) , as.matrix(Zi%*%t(Zi))))
            
            
            sigZi <- kronecker(cholinv(prior_sigma), Zi)
            errvec <- matrix(as.vector(t(Yi - prior_beta%*%Xi)), byrow = T, ncol = 1)
            psi_hat <- Ui%*%sigZi%*%errvec
            
            post_psivec <- rmvn(1, as.vector(psi_hat), Ui)
            post_psi[,k] <- post_psivec
          }
          prior_psi_star <- post_psi
          post_psi <- matrix(as.vector(post_psi), nrow = d)
          psitrack[i+1,] <- as.vector(post_psi)
          prior_psi <- post_psi
        } else{
          #To troubleshoot: keep Psi fixed
          prior_psi_star <- matrix(as.vector(prior_psi), ncol = m) #Not entirely sure this will work with multiple random effects
          post_psi <- prior_psi
          psitrack[i+1,] <- as.vector(post_psi)
        }
        
      }
      
      #Phi|Psi
      
      if(fixPhi == FALSE){
        Lambda_2n <- Lambda_2 + tcrossprod(prior_psi_star)
        post_phi <- riwish(v = v_2n, S = Lambda_2n)
        phitrack[i+1,] <- as.vector(post_phi)
        prior_phi <- post_phi
      } else{
        #To troubleshoot: keep Phi fixed
        post_phi <- prior_phi
        phitrack[i+1,] <- as.vector(prior_phi)
      }
      
      if(fixBeta == FALSE){
        #Y given X sums of squares matrix for Sigma constructed using beta_hat and updated Psi
        if(block == FALSE){
          S_ygx <- tcrossprod((Y - tcrossprod(beta_hat, t(X)) - prior_psi%*%as.matrix(t(Z))))
        } else{
          S_ygx <- tcrossprod((Y - tcrossprod(prior_beta, t(X)) - prior_psi%*%as.matrix(t(Z))))
        }
      } else{
        S_ygx <- tcrossprod((Y - tcrossprod(prior_beta, t(X)) - prior_psi%*%as.matrix(t(Z))))
      }
      
      
      #Sigma| Y, X, Beta, Psi
      
      if(fixSigma == FALSE){
        Lambda_1n <- Lambda_1 + S_ygx
        post_sigma <- riwish(v = v_1n, S = Lambda_1n)
        sigmatrack[i+1,] <- as.vector(post_sigma)
        prior_sigma <- post_sigma
      } else{
        #To Troubleshoot: keep Sigma fixed
        post_sigma <- prior_sigma
        sigmatrack[i+1,] <- as.vector(post_sigma)
      }
      
      #Beta| Y, X, Z, Psi, Sigma
      
      if(block == FALSE){
        
        if(fixBeta == FALSE){
          post_beta_vec_cov <- kronecker(cholinv(S_xx), prior_sigma)
          post_beta_0 <- beta_hat
          post_beta_0_vec <- as.vector(post_beta_0)
          post_beta_vec <- rmvn(n = 1, mu = post_beta_0_vec, sigma = post_beta_vec_cov)
          post_beta <- matrix(post_beta_vec, nrow = d) #should be d by p
          betatrack[i+1,] <- as.vector(post_beta)
          prior_beta <- post_beta
        } else{
          #To Troubleshoot: keep Beta fixed
          post_beta <- prior_beta
          betatrack[i+1,] <- as.vector(post_beta)
        }
        
      }
      
      if(fixBeta == FALSE){
        rm(S_yx, margS_yx, beta_hat, margbeta_hat) #clean up memory
      }
      
      #Update log-odds means, Mu
      prior_mu <- tcrossprod(prior_beta, t(X)) + prior_psi%*%as.matrix(t(Z))
      
    }
    
    paths <- list(Y = ytrack, Sigma = sigmatrack, Beta = betatrack, Psi = psitrack, Phi = phitrack, MHRatio = MHacceptancelist)
    return(paths)
    
  }
  
  #Fixed Model
  if(mixed == FALSE){
    
    #Some dimensions
    p <- ncol(X) #number of covariates, including intercept
    d <- ncol(W) - 1 #number of multivariate categories
    N <- nrow(X) #number of total observations
    
    #Algorithm uses a wide X
    X <- t(X) #should be p by N
    
    #Get initial estimates of the log-odds
    if(fixY == FALSE){
      Y <- t(newinitialY(W, base)) #should be d by N
    } else{
      Y <- t(trY)
    }
    
    #Beta set-up
    if(fixBeta == FALSE){
      prior_beta <- t(tcrossprod(cholinv(tcrossprod(X)), tcrossprod(Y, X)))
      betatrack <- matrix(0, nrow = iter + 1, ncol = d*p)
      betatrack[1,] <- as.vector(prior_beta)
      if(startru == TRUE){
        prior_beta <- trBeta
        betatrack <- matrix(0, nrow = iter + 1, ncol = d*p)
        betatrack[1,] <- as.vector(prior_beta)
      }
    } else{
      prior_beta <- trBeta
      betatrack <- matrix(0, nrow = iter + 1, ncol = d*p)
      betatrack[1,] <- as.vector(prior_beta)
    }
    
    #Sigma set-up
    #Hyperparameters for prior distribution of Sigma
    v_1 <- d + 1
    Lambda_1 <- diag(d)
    
    if(fixSigma == FALSE){
      prior_sigma <- diag(d)
      sigmatrack <- matrix(0, nrow = iter + 1, ncol = length(prior_sigma))
      sigmatrack[1,] <- as.vector(prior_sigma)
      if(startru == TRUE){
        prior_sigma <- trSigma
        sigmatrack <- matrix(0, nrow = iter + 1, ncol = length(prior_sigma))
        sigmatrack[1,] <- as.vector(prior_sigma)
      }
    } else{
      prior_sigma <- trSigma
      sigmatrack <- matrix(0, nrow = iter + 1, ncol = length(prior_sigma))
      sigmatrack[1,] <- as.vector(prior_sigma)
    }
    
    #Get initial log-odds means, and then set starting values of Y to that (except doing that leads to horrible MH acceptance ratio when moving to more categories)
    prior_mu <- tcrossprod(prior_beta, t(X))
    
    #Starting values for Y
    if(fixY == FALSE){
      if(whichYstart == "mean"){
        Y <- prior_mu
      }
    } else{
      Y <- t(trY)
    }
    
    ytrack <- matrix(0, nrow = iter + 1, ncol = d*N)
    ytrack[1,] <- as.vector(Y)
    
    #Certain values will stay fixed throughout the Gibbs chain
    v_1n <- v_1 + N #posterior df of Sigma
    S_xx <- tcrossprod(X) #X sums of squares matrix
    
    #Finally, also want to keep track of MH acceptance probabilities for proposing Y
    MHacceptancelist <- matrix(0, nrow = iter, ncol = mhiter)
    
    #Full posterior log-likelihood and error term tracking can be recovered from the resulting output matrices, and to save computation time are ommitted from the main function
    
    #The MH-within-Gibbs Algorithm
    for(i in 1:iter){
      
      #Metropolis-Hastings step
      #Y|W, X, Z, Beta, Psi, Sigma, Phi
      
      if(fixY == FALSE){
        for(j in 1:mhiter){
          Y <- t(Y) #put it long ways for proposing
          Mu <- t(prior_mu) #put it long ways for proposing
          wmu <- cbind(W, Mu) #for the proposal function
          wmuy <- cbind(wmu, Y)
          
          if(proposal == "norm"){
            Y_new <- Y + rmvn(N, integer(length = d), mhscale*prior_sigma) #propose centered at old Y
          }
          
          if(proposal == "beta"){
            Pstar <- t(apply(X = Y, MARGIN = 1, FUN = ytopstar)) #old Pstar
            Pstar_new <- t(apply(X = wmu, MARGIN = 1, FUN = betapropdist, Sigma = mhscale*prior_sigma))
            Y_new <- t(apply(X = Pstar_new, MARGIN = 1, FUN = pstartoy)) #newY based on proposed probabilities
            wmuy_new <- cbind(wmu, Y_new)
            wmuPstar <- cbind(wmu, Pstar)
            wmuPstar_new <- cbind(wmu, Pstar_new)
          }
          
          if(proposal == "normbeta"){
            Y_new <- t(apply(X = wmu, MARGIN = 1, FUN = normbetapropdist, Sigma = mhscale*prior_sigma)) #Fixed Normal proposal based on logit(Beta) expected values and variances
            wmuy_new <- cbind(wmu, Y_new)
          }
          
          sexpY <- rowSums(exp(Y))
          sexpYn <- rowSums(exp(Y_new))
          
          newnormpart <- .5*tcrossprod(tcrossprod((Y_new - Mu) , t(cholinv(prior_sigma))), (Y_new - Mu))
          oldnormpart <- .5*tcrossprod(tcrossprod((Y - Mu) , t(cholinv(prior_sigma))), (Y - Mu))
          newloglike <- rowSums(as.matrix(W[,1:d]*Y_new[,1:d])) - rowSums(W*(log1p(sexpYn))) - newnormpart[cbind(1:nrow(newnormpart), 1:nrow(newnormpart))]
          oldloglike <- rowSums(as.matrix(W[,1:d]*Y[,1:d])) - rowSums(W*(log1p(sexpY))) - oldnormpart[cbind(1:nrow(oldnormpart), 1:nrow(oldnormpart))]
          
          rm(newnormpart, oldnormpart) #remove from memory to help efficiency
          
          if(proposal == "norm"){
            #Random walk is symmetric, Metropolis, not MH
            logptolike <- 0
            logotplike <- 0
            
            rm(wmu, wmuy)
          }
          
          if(proposal == "beta"){
            logptolike <- apply(X = wmuPstar, MARGIN = 1, FUN = betaloglike, Sigma = mhscale*prior_sigma) #going to old p_star under proposal distribution, including Jacobian
            logotplike <- apply(X = wmuPstar_new, MARGIN = 1, FUN = betaloglike, Sigma = mhscale*prior_sigma) #going to proposed p_star under proposal, including Jacobian
            
            rm(wmu, wmuy, wmuy_new, wmuPstar, wmuPstar_new)
          }
          
          if(proposal == "normbeta"){
            logptolike <- apply(X = wmuy, MARGIN = 1, FUN = normbetaloglike, Sigma = mhscale*prior_sigma)
            logotplike <- apply(X = wmuy_new, MARGIN = 1, FUN = normbetaloglike, Sigma = mhscale*prior_sigma)
            
            rm(wmu, wmuy, wmuy_new)
          }
          
          ratio <- newloglike - oldloglike + logptolike - logotplike
          temp <- runif(N)
          difr <- as.numeric((ratio - log(temp)) >= 0)
          MHacceptancelist[i,j] <- mean(difr)
          Y <- difr*Y_new + (1 - difr)*Y #accept or reject samples to get draw
          Y <- t(Y) #put it back for the rest of the Gibbs
        }
        
        ytrack[i+1,] <- as.vector(Y)
      } else{
        #To Troubleshoot: keep Y fixed
        ytrack[i+1,] <- as.vector(Y)
      }
      
      #Using sampled Y, create necessary sums of squares matrices and OLS mean estimates for:
      #individual update of Beta
      if(fixBeta == FALSE){
        S_yx <- tcrossprod(Y, X)
        beta_hat <- tcrossprod(S_yx, t(cholinv(S_xx)))
      }
      
      #If Fixed, there should be no block update
      
      if(fixBeta == FALSE){
        #Y given X sums of squares matrix for Sigma constructed using beta_hat and updated Psi
        S_ygx <- tcrossprod(Y - tcrossprod(beta_hat, t(X)))
      } else{
        S_ygx <- tcrossprod(Y - tcrossprod(prior_beta, t(X)))
      }
      
      #Sigma| Y, X, Beta, Psi
      
      if(fixSigma == FALSE){
        Lambda_1n <- cholinv(Lambda_1) + S_ygx
        post_sigma <- riwish(v = v_1n, S = Lambda_1n)
        sigmatrack[i+1,] <- as.vector(post_sigma)
        prior_sigma <- post_sigma
      } else{
        #To Troubleshoot: keep Sigma fixed
        post_sigma <- prior_sigma
        sigmatrack[i+1,] <- as.vector(post_sigma)
      }
      
      #Beta| Y, X, Z, Psi, Sigma
      
      if(fixBeta == FALSE){
        post_beta_vec_cov <- kronecker(cholinv(S_xx), prior_sigma)
        post_beta_0 <- beta_hat
        post_beta_0_vec <- as.vector(post_beta_0)
        post_beta_vec <- rmvn(n = 1, mu = post_beta_0_vec, sigma = post_beta_vec_cov)
        post_beta <- matrix(post_beta_vec, nrow = d) #should be d by p
        betatrack[i+1,] <- as.vector(post_beta)
        prior_beta <- post_beta
      } else{
        #To Troubleshoot: keep Beta fixed
        post_beta <- prior_beta
        betatrack[i+1,] <- as.vector(post_beta)
      }
      
      if(fixBeta == FALSE){
        rm(S_yx, beta_hat) #clean up memory
      }
      
      #Update log-odds means, Mu
      prior_mu <- tcrossprod(prior_beta, t(X))
      
    }
    
    paths <- list(Y = ytrack, Sigma = sigmatrack, Beta = betatrack, MHRatio = MHacceptancelist)
    return(paths)
    
  }
  
  
}

#Once the algorithm is completed, can input the output matrices into this function to look at the convergence of the full posterior log-likelihood
logliketrack <- function(W, X, Z, randeff, ytrack, betatrack, psitrack, sigmatrack, phitrack, block = TRUE, test = FALSE){
  
  
  #Some dimensions
  p <- ncol(X) #number of covariates, including intercept
  d <- ncol(W) - 1 #number of multivariate categories
  N <- nrow(X) #number of total observations
  q <- randeff #number of random player effects (default is simply an intercept)
  m <- ncol(Z)/randeff #number of players
  
  #Algorithm uses a wide X
  X <- t(X) #should be p by N
  
  #Hyperparameters for prior distribution of Phi
  v_2 <- d*q
  Lambda_2 <- diag(d*q) + matrix(1, nrow = d*q)%*%matrix(1, ncol = d*q)
  
  #Hyperparameters for prior distribution of Sigma
  v_1 <- d
  Lambda_1 <- diag(d) + matrix(1, nrow = d)%*%matrix(1, ncol = d)
  
  full_post_log_like_track <- matrix(0, nrow = nrow(ytrack) - 1, ncol = 6)
  
  if(test == TRUE){
    browser()
  }
  
  getdatlike <- function(wpy_vec){
    sexpy <- 1 + sum(exp(wpy_vec[(d+2):length(wpy_vec)]))
    yprob <- c((exp(wpy_vec[(d+2):length(wpy_vec)])/sexpy), (1/sexpy))
    result <- dmultinom(wpy_vec[1:(d+1)], prob = yprob, log = T)
    return(result)
  }
  
  getdatlike2 <- function(psii, phii){
    result <- dmvn(psii, mu = rep(0, d), sigma = phii, log = TRUE)
    return(result)
  }
  
  for(i in 2:nrow(ytrack)){
    
    Y <- matrix(ytrack[i,], ncol = N)
    WpY <- rbind(t(W), Y)
    
    #Multinomial part: W|Y
    #      ll_wgy <- 0
    #      for(l in 1:nrow(W)){
    #        sumexpyip1 <- 1 + sum(exp(Y[,l]))
    #        y_prob <- unlist(c((exp(Y[,l])/sumexpyip1), (1/sumexpyip1)))
    #        ll_wgyi <- dmultinom(W[l,], prob = y_prob, log = TRUE)
    #        ll_wgy <- ll_wgy + ll_wgyi
    #      }
    
    #Re-Code to be more efficient
    ll_wgy <- sum(apply(WpY, 2, getdatlike))
    
    ##Y's normal part: Y|Psi, Beta, Sigma, X
    pmi <- tcrossprod(matrix(betatrack[i-1,], byrow = F, nrow = d), t(X)) + matrix(psitrack[i-1,], byrow = F, nrow = d)%*%as.matrix(t(Z))
    ll_ygmat <- -.5*log(det(2*pi*matrix(sigmatrack[i-1,], byrow = F, ncol = d))) - .5*tcrossprod(crossprod((Y - pmi), cholinv(matrix(sigmatrack[i-1,], byrow = F, ncol = d))), t(Y - pmi))
    ll_ygpbs <- sum(ll_ygmat[cbind(1:nrow(ll_ygmat), 1:nrow(ll_ygmat))])
    
    ##Psi's normal part: Psi|Phi
    #ll_pgp <- dmvn(psitrack[i,], mu = rep(0, d*q*m), sigma = kronecker(diag(m), matrix(phitrack[i-1,], byrow = F, ncol = d*q)), log = TRUE)
    psii <- matrix(psitrack[i,], ncol = d, byrow = T)
    
    #Re-Code to be more efficient
    ll_pgp <- sum(apply(psii, 1, getdatlike2, phii = matrix(phitrack[i-1,], byrow = T, ncol = d)))
    
    ##Beta's part: Beta|Sigma (depends on block)
    #Block update Beta prior is uniform
    if(block == TRUE){
      ll_bgs <- log(1)
    } else{
      #Individual update, Beta prior is normal
      ll_bgs <- dmvn(betatrack[i,], mu = rep(0, d*p), sigma = kronecker(cholinv(tcrossprod(X)), matrix(sigmatrack[i-1,], byrow = F, ncol = d)), log = TRUE)
    }
    
    ##Phi's inverse-Wishart part: Phi
    ll_p <- log(diwish(matrix(phitrack[i,], ncol = d*q), v = v_2, S = Lambda_2))
    
    ##Sigma's inverse-wishart part: Sigma
    ll_s <- log(diwish(matrix(sigmatrack[i,], byrow = F, ncol = d), v = v_1, S = Lambda_1))
    
    ##Put them all together and what do you get? The full posterior log likelihood
    full_post_log_like_track[i-1,] <- c(ll_wgy, ll_ygpbs, ll_pgp, ll_bgs, ll_p, ll_s)
    
  }
  
  return(full_post_log_like_track)
  
}

#Trueloglikelihood
truloglike <- function(W, X, Z, randeff, y, beta, psi, sigma, phi, block = TRUE, test = FALSE){
  
  
  #Some dimensions
  p <- ncol(X) #number of covariates, including intercept
  d <- ncol(W) - 1 #number of multivariate categories
  N <- nrow(X) #number of total observations
  q <- randeff #number of random player effects (default is simply an intercept)
  m <- ncol(Z)/randeff #number of players
  
  #Algorithm uses a wide X
  X <- t(X) #should be p by N
  
  #Hyperparameters for prior distribution of Phi
  #v_2 <- d*q
  #Lambda_2 <- diag(d*q) + matrix(1, nrow = d*q)%*%matrix(1, ncol = d*q)
  #Under Gelman
  v_2 <- d*q + 2
  Lambda_2 <- diag(d*q)*(1/(2*.001))
  
  #Hyperparameters for prior distribution of Sigma
  #v_1 <- d
  #Lambda_1 <- diag(d) + matrix(1, nrow = d)%*%matrix(1, ncol = d)
  #Under Gelman
  v_1 <- d + 2
  Lambda_1 <- diag(d)*(1/(2*.001))
  
  #Multinomial part: W|Y
  ll_wgy <- 0
  for(l in 1:nrow(W)){
    sumexpyip1 <- 1 + sum(exp(y[l,]))
    y_prob <- unlist(c((exp(y[l,])/sumexpyip1), (1/sumexpyip1)))
    ll_wgyi <- dmultinom(W[l,], prob = y_prob, log = TRUE)
    ll_wgy <- ll_wgy + ll_wgyi
  }
  
  if(test == TRUE){
    browser()
  }
  
  ##Y's normal part: Y|Psi, Beta, Sigma, X
  pmi <- tcrossprod(t(beta), t(X)) + t(psi)%*%as.matrix(t(Z))
  ll_ygmat <- -.5*log(det(2*pi*sigma)) - .5*tcrossprod(crossprod((t(y) - pmi), cholinv(sigma)), t(t(y) - pmi))
  ll_ygpbs <- sum(ll_ygmat[cbind(1:nrow(ll_ygmat), 1:nrow(ll_ygmat))])
  
  ##Psi's normal part: Psi|Phi
  ll_pgp <- dmvn(as.vector(t(psi)), mu = rep(0, d*q*m), sigma = kronecker(diag(m), matrix(as.vector(phi), byrow = F, ncol = d*q)), log = TRUE)
  
  ##Beta's part: Beta|Sigma (depends on block)
  #Block update Beta prior is uniform
  if(block == TRUE){
    ll_bgs <- log(1)
  } else{
    #Individual update, Beta prior is normal
    ll_bgs <- dmvn(as.vector(t(beta)), mu = rep(0, d*p), sigma = kronecker(cholinv(tcrossprod(X)), matrix(as.vector(sigma), byrow = F, ncol = d)), log = TRUE)
  }
  
  ##Phi's inverse-Wishart part: Phi
  ll_p <- log(diwish(matrix(as.vector(phi), ncol = d*q), v = v_2, S = Lambda_2))
  
  ##Sigma's inverse-wishart part: Sigma
  ll_s <- log(diwish(matrix(as.vector(sigma), byrow = F, ncol = d), v = v_1, S = Lambda_1))
  
  
  return(c(ll_wgy, ll_ygpbs, ll_pgp, ll_bgs, ll_p, ll_s))
  
}

#Mixed Effects Multinomial Logit Gibbs Sampler, with proposal covariance matrix based on Polya-Gamma auxiliary variables
multmixgibbs_b0prior <- function(W, X, Z, randeff = 1, base, iter = 1000, mhiter_psi = 1, mhiter_beta = 1, mhscale_psi = 1, mhscale_int = 1, mhscale_beta = 1, r.seed = 42, trY = NULL, trBeta = NULL, trPsi = NULL, trPhi = NULL, startru = FALSE, startBetatru = FALSE, startPsitru = FALSE, fixIntercept = FALSE, fixBeta = FALSE, fixPsi = FALSE, fixPhi = FALSE, whichYstart = c("mean", "inity"), test = FALSE){
  
  #Set seed for replicability
  set.seed(r.seed)
  
  #Some dimensions
  p <- ncol(X) #number of covariates, including intercept
  d <- ncol(W) - 1 #number of multivariate categories
  N <- nrow(X) #number of total observations
  q <- randeff #number of random player effects (default is simply an intercept)
  m <- ncol(Z)/randeff #number of players
  
  pa <- rowSums(W)
  
  #Algorithm uses a wide X
  X <- t(X) #should be p by N
  
  #Get initial estimates of the log-odds
  Y <- t(newinitialY(W, base)) #should be d by N
  if(startru == TRUE){
    Y <- t(trY)
  }
  
  #Beta set-up
  if(fixBeta == FALSE){
    initial_beta <- t(tcrossprod(cholinv(tcrossprod(X)), tcrossprod(Y, X)))
    prior_beta <- initial_beta[,-1]
    betatrack <- matrix(0, nrow = iter + 1, ncol = d*(p-1))
    betatrack[1,] <- as.vector(prior_beta)
    if(startBetatru == TRUE){
      prior_beta <- trBeta[,-1]
      betatrack <- matrix(0, nrow = iter + 1, ncol = d*(p-1))
      betatrack[1,] <- as.vector(prior_beta)
    }
  } else{
    prior_beta <- trBeta[,-1]
    betatrack <- matrix(0, nrow = iter + 1, ncol = d*(p-1))
    betatrack[1,] <- as.vector(prior_beta)
  }
  
  if(fixIntercept == FALSE){
    initial_beta <- t(tcrossprod(cholinv(tcrossprod(X)), tcrossprod(Y, X)))
    prior_intercept <- initial_beta[,1]
    intercepttrack <- matrix(0, nrow = iter + 1, ncol = d)
    intercepttrack[1,] <- as.vector(prior_intercept)
    if(startru == TRUE){
      prior_intercept <- trBeta[,1]
      intercepttrack <- matrix(0, nrow = iter + 1, ncol = d)
      intercepttrack[1,] <- as.vector(prior_intercept)
    }
  } else{
    prior_intercept <- trBeta[,1]
    intercepttrack <- matrix(0, nrow = iter + 1, ncol = d)
    intercepttrack[1,] <- prior_intercept
  }
  
  #Get initial log-odds means, and then set starting values of Y to that (except doing that leads to horrible MH acceptance ratio when moving to more categories)
  prior_mu <- tcrossprod(prior_beta, t(X[-1,]))
  
  #Auxiliary variables for the Polya-Gamma framework, call them Gamma, there are N*d of them
  #Posterior given Beta is gamma_ijk ~ PG(pa_ij = sum(w_ijk), nu_ijk = y_ijk - c_ijk), where c_ijk = log(sum (notk) exp(y_ijk)); just sample from that for initial values
  
  #Set up C_ijk
  c_ijk <- matrix(0, ncol = d, nrow = N)
  for(k in 1:d){
    if(d == 2){
      c_ijk[,k] <- log(1 + exp(prior_mu[-k,]))
    }else{
      c_ijk[,k] <- log(1 + colSums(exp(prior_mu[-k,])))
    }
  }
  
  #For some reason BayesLogit (or anything that allows me to sample a polya-gamma random variate) is not available for R on my laptop, so we'll use a different proposal covariance structure for the time being (it shouldn't affect the likelihood)
  
  #Auxiliary variables
  prior_gamma_ijk <- matrix(0, nrow = d, ncol = N)
  for(k in 1:d){
    if(d == 2){
      prior_gamma_ijk[k,] <- apply(rbind((prior_mu[k,] - c_ijk[,k]), pa), 2, function(nupa) { rpg(num = 1, h = nupa[2], z = nupa[1]) })
    }else{
      prior_gamma_ijk[k,] <- apply(rbind((prior_mu[k,] - c_ijk[,k]), pa), 2, function(nupa) { rpg(num = 1, h = nupa[2], z = nupa[1]) })
    }
  }
  
  #Psi set-up
  if(fixPsi == FALSE){
    prior_psi <- matrix(prior_intercept, nrow = d, ncol = q*m, byrow = F)
    psitrack <- matrix(0, nrow = iter + 1, ncol = length(prior_psi))
    psitrack[1,] <- as.vector(prior_psi)
    if(startPsitru == TRUE){
      prior_psi <- t(trPsi) #should be d by q*m
      psitrack <- matrix(0, nrow = iter + 1, ncol = length(prior_psi))
      psitrack[1,] <- as.vector(prior_psi)
    }
  } else{
    prior_psi <- t(trPsi)
    psitrack <- matrix(0, nrow = iter + 1, ncol = length(prior_psi))
    psitrack[1,] <- as.vector(prior_psi)
  }
  
  #Starting values for Y
  if(whichYstart == "mean"){
    Y <- prior_mu + prior_psi%*%as.matrix(t(Z))
    if(startru == TRUE){
      Y <- t(trY)
    }
  }
  
  #Phi set-up
  #Hyperparameters for prior distribution of Phi
  v_2 <- d*q
  Lambda_2 <- diag(d*q) + matrix(1, nrow = d*q)%*%matrix(1, ncol = d*q)
  
  if(fixPhi == FALSE){
    prior_phi <- diag(d*q)
    phitrack <- matrix(0, nrow = iter + 1, ncol = length(prior_phi))
    phitrack[1,] <- as.vector(prior_phi)
    if(startru == TRUE){
      prior_phi <- trPhi
      phitrack <- matrix(0, nrow = iter + 1, ncol = length(prior_phi))
      phitrack[1,] <- as.vector(prior_phi)
    }
  } else{
    prior_phi <- trPhi
    phitrack <- matrix(0, nrow = iter + 1, ncol = length(prior_phi))
    phitrack[1,] <- as.vector(prior_phi)
  }
  
  #Am going to keep track of the Y's only after the update Psi step (can always go back and look at what they are after the Beta step by adding XBeta + Psi for an iteration)
  ytrack <- matrix(0, nrow = iter + 1, ncol = d*N)
  ytrack[1,] <- as.vector(Y)
  
  #Certain values will stay fixed throughout the Gibbs chain
  v_2n <- v_2 + m #posterior df of Phi
  
  #Finally, also want to keep track of MH acceptance probabilities for proposing Y
  MHacceptancelist <- matrix(0, nrow = iter, ncol = mhiter_psi)
  
  MHacceptancelist_beta <- matrix(0, nrow = iter, ncol = d)
  
  #Full posterior log-likelihood and error term tracking can be recovered from the resulting output matrices, and to save computation time are ommitted from the main function
  
  #The MH-within-Gibbs Algorithm
  for(i in 1:iter){
    
    #This if statement can be moved around for troubleshooting at different points
    if(test == TRUE){
      browser()
    }
    
    #        if(i %in% seq(100, iter, round(iter/100))){
    #          print(i)
    #        }
    
    #Metropolis-Hastings step
    #Y|W, X, Z, Beta, Phi
    
    if(fixPsi == FALSE){
      for(j in 1:mhiter_psi){
        Y <- t(Y) #put it long ways for proposing
        Mu <- t(prior_mu) #put it long ways for proposing
        playerid <- apply(Z, 1, function(x) which(x == 1))
        seasons <- as.vector(table(playerid)) #will only work for single random intercept
        fullPsi <- Y - Mu
        Psi_mean <- matrix(prior_intercept, ncol = d, byrow = T, nrow = m)
        oldPsi <- as.matrix(t(Z)%*%fullPsi/seasons)
        
        #I think, in this situation, only proposing from a normal is possible...
        
        psi_new <- oldPsi + rmvn(m, integer(length = d), mhscale_psi*prior_phi) #propose centered at old psi
        
        makenewfullPsi <- function(i, psi){
          rep(psi[i,], seasons[i])
        }
        
        newfullPsi <- matrix(unlist(sapply(1:m, makenewfullPsi, psi = psi_new)), ncol = d, byrow = T)
        
        Y_new <- Mu + newfullPsi
        
        sexpY <- rowSums(exp(Y))
        sexpYn <- rowSums(exp(Y_new))
        
        newnormpart <- .5*tcrossprod(tcrossprod((psi_new - Psi_mean) , t(cholinv(prior_phi))), (psi_new - Psi_mean))
        oldnormpart <- .5*tcrossprod(tcrossprod((oldPsi - Psi_mean) , t(cholinv(prior_phi))), (oldPsi - Psi_mean))
        newloglike <- rowSums(as.matrix(W[,1:d]*Y_new[,1:d])) - rowSums(W*(log1p(sexpYn)))
        oldloglike <- rowSums(as.matrix(W[,1:d]*Y[,1:d])) - rowSums(W*(log1p(sexpY)))
        
        #No longer can I accept/reject player-seasons, I only accept/reject players....
        newloglike_player <- unlist(lapply(split(newloglike, playerid), sum)) - newnormpart[cbind(1:nrow(newnormpart), 1:nrow(newnormpart))]
        oldloglike_player <- unlist(lapply(split(oldloglike, playerid), sum)) - oldnormpart[cbind(1:nrow(oldnormpart), 1:nrow(oldnormpart))]
        
        ratio_player <- newloglike_player - oldloglike_player
        
        temp_player <- runif(m)
        difr_player <- as.numeric((ratio_player - log(temp_player)) >= 0)
        MHacceptancelist[i,j] <- mean(difr_player)
        post_psi <- difr_player*psi_new + (1 - difr_player)*oldPsi
        prior_psi <- t(post_psi)
        fullPsi <- matrix(unlist(sapply(1:m, makenewfullPsi, psi = post_psi)), ncol = d, byrow = T)
        
        Y <- Mu + fullPsi
        Y <- t(Y) #put it back for the rest of the Gibbs
        
      }
      
      psitrack[i+1,] <- as.vector(prior_psi)
      ytrack[i+1,] <- as.vector(Y)
    } else{
      #To Troubleshoot: keep Psi fixed 
      playerid <- apply(Z, 1, function(x) which(x == 1))
      seasons <- as.vector(table(playerid)) #will only work for single random intercept
      makenewfullPsi <- function(i, psi){
        rep(psi[i,], seasons[i])
      }
      fullPsi <- matrix(unlist(sapply(1:m, makenewfullPsi, psi = trPsi)), ncol = d, byrow = T)
      psitrack[i+1,] <- as.vector(prior_psi)
      ytrack[i+1,] <- as.vector(Y)
    }
    
    #Phi|Psi
    
    if(fixPhi == FALSE){
      prior_psi_star <- prior_psi
      Psi_mean <- matrix(rowMeans(prior_psi_star), nrow = d, byrow = F, ncol = m)
      Lambda_2n <- Lambda_2 + tcrossprod(prior_psi_star - Psi_mean)
      post_phi <- riwish(v = v_2n, S = Lambda_2n)
      phitrack[i+1,] <- as.vector(post_phi)
      prior_phi <- post_phi
    } else{
      #To troubleshoot: keep Phi fixed
      post_phi <- prior_phi
      phitrack[i+1,] <- as.vector(prior_phi)
    }
    
    #Intercept|Psi, Phi
    if(fixIntercept == FALSE){
      post_intercept <- rmvn(1, rowMeans(prior_psi), prior_phi/m)
      intercepttrack[i+1,] <- as.vector(post_intercept)
      prior_intercept <- as.vector(post_intercept)
    } else{
      post_intercept <- prior_intercept
      intercepttrack[i+1,] <- as.vector(prior_intercept)
    }
    
    
    #Need to write an entirely new Beta step given data (since no longer conditioning on Y), which will likely be a Metropolis step; akin to the update of the Betas in regular multinomial logit
    #To perhaps speed up mixing, let's do an independent update of each Beta vector based on the conditional distribution in the Polya-Gamma paper?
    
    if(fixBeta == FALSE){
      Y <- t(Y) #put it long ways for proposing
      Mu <- t(prior_mu) #put it long ways for proposing
      oldBeta <- prior_beta #should be the correct dimension (Beta is d*(p-1), since we do the intercept separately)
      
      for(k in 1:d){
        kappa_k <- W[,k] - rowSums(W)/2
        omega_k <- diag(prior_gamma_ijk[k,]) #N by N matrix
        post_beta_cova_noint <- cholinv(X[-1,]%*%omega_k%*%t(X[-1,])) #p-1 by p-1 matrix
        
        beta_new_k <- oldBeta[k,] + rmvn(1, integer(length = p - 1), mhscale_beta*post_beta_cova_noint)
        
        oldnu <- Y[,k] - c_ijk[,k]
        newnu <- t(X[-1,])%*%matrix(beta_new_k, nrow = p-1) + fullPsi[,k] - c_ijk[,k]
        
        #Symmetric proposal means Metropolis
        newloglike_beta <- sum(W[,k]*newnu - W[,k]*log1p(exp(newnu)) - rowSums(W[,-k])*log1p(exp(newnu)))
        oldloglike_beta <- sum(W[,k]*oldnu - W[,k]*log1p(exp(oldnu)) - rowSums(W[,-k])*log1p(exp(oldnu)))
        
        ratio_beta <- newloglike_beta - oldloglike_beta
        
        temp_beta <- runif(1)
        difr_beta <- as.numeric((ratio_beta - log(temp_beta)) >= 0)
        MHacceptancelist_beta[i,k] <- difr_beta
        
        post_beta_k <- difr_beta*beta_new_k + (1 - difr_beta)*oldBeta[k,]
        prior_beta[k,] <- post_beta_k
        
      }
      
      betatrack[i+1,] <- as.vector(prior_beta)
      
      Y <- t(tcrossprod(prior_beta, t(X[-1,]))) + fullPsi
      Y <- t(Y) #put it back for the rest of the Gibbs
      
      #Update C_ijk and Auxiliary variables
      c_ijk <- matrix(0, ncol = d, nrow = N)
      for(k in 1:d){
        if(d == 2){
          c_ijk[,k] <- log(1 + exp(Y[-k,]))
        }else{
          c_ijk[,k] <- log(1 + colSums(exp(Y[-k,])))
        }
      }
      
      #Auxiliary variables
      prior_gamma_ijk <- matrix(0, nrow = d, ncol = N)
      for(k in 1:d){
        if(d == 2){
          prior_gamma_ijk[k,] <- apply(rbind((Y[k,] - c_ijk[,k]), pa), 2, function(nupa) { rpg(num = 1, h = nupa[2], z = nupa[1]) })
        }else{
          prior_gamma_ijk[k,] <- apply(rbind((Y[k,] - c_ijk[,k]), pa), 2, function(nupa) { rpg(num = 1, h = nupa[2], z = nupa[1]) })
        }
      }
      
    } else{
      #To Troubleshoot: keep Beta fixed
      betatrack[i+1,] <- as.vector(prior_beta)
    }
    
    #Update log-odds means, Mu
    prior_mu <- tcrossprod(prior_beta, t(X[-1,]))
    
  }
  
  paths <- list(Y = ytrack, Intercept = intercepttrack, Beta = betatrack, Psi = psitrack, Phi = phitrack, MHRatio_Psi = MHacceptancelist, MHRatio_Beta = MHacceptancelist_beta)
  return(paths)
  
}

#Need to check the convergence of the complete data loglikelihood (M-Logit)
multlogit_logliketrack <- function(W, X, Z, randeff, ytrack, mutrack, betatrack, psitrack, phitrack, test = FALSE){
  
  
  #Some dimensions
  p <- ncol(X) #number of covariates, including intercept
  d <- ncol(W) - 1 #number of multivariate categories
  N <- nrow(X) #number of total observations
  q <- randeff #number of random player effects (default is simply an intercept)
  m <- ncol(Z)/randeff #number of players
  
  #Algorithm uses a wide X
  X <- t(X) #should be p by N
  
  #Hyperparameters for prior distribution of Phi
  v_2 <- d*q
  Lambda_2 <- diag(d*q) + matrix(1, nrow = d*q)%*%matrix(1, ncol = d*q)
  
  full_post_log_like_track <- matrix(0, nrow = nrow(betatrack) - 1, ncol = 5)
  
  getdatlike <- function(wpy_vec){
    sexpy <- 1 + sum(exp(wpy_vec[(d+2):length(wpy_vec)]))
    yprob <- c((exp(wpy_vec[(d+2):length(wpy_vec)])/sexpy), (1/sexpy))
    result <- dmultinom(wpy_vec[1:(d+1)], prob = yprob, log = T)
    return(result)
  }
  
  getdatlike2 <- function(psii, mui, phii){
    result <- dmvn(psii, mu = mui, sigma = phii, log = TRUE)
    return(result)
  }
  
  for(i in 2:nrow(betatrack)){
    
    Y <- matrix(ytrack[i,], ncol = N)
    WpY <- rbind(t(W), Y)
    
    #Multinomial part: W|Y
    #ll_wgy <- 0
    #for(l in 1:nrow(W)){
    #  sumexpyip1 <- 1 + sum(exp(Y[,l]))
    #  y_prob <- unlist(c((exp(Y[,l])/sumexpyip1), (1/sumexpyip1)))
    #  ll_wgyi <- dmultinom(W[l,], prob = y_prob, log = TRUE)
    #  ll_wgy <- ll_wgy + ll_wgyi
    #}
    
    #Re-Code for Efficiency
    ll_wgy <- sum(apply(WpY, 2, getdatlike))
    
    if(test == TRUE){
      browser()
    }
    
    ##Psi's normal part: Psi|Mu, Phi
    #ll_pgp <- dmvn(psitrack[i,], mu = rep(mutrack[i-1,], m), sigma = kronecker(diag(m), matrix(phitrack[i-1,], byrow = F, ncol = d*q)), log = TRUE)
    psii <- matrix(psitrack[i,], ncol = d, byrow = T)
    
    #Re-Code for Efficiency
    ll_pgp <- sum(apply(psii, 1, getdatlike2, mui = mutrack[i-1,], phii = matrix(phitrack[i-1,], byrow = T, ncol = d)))
    
    ##Mu's part: uniform
    ll_m <- log(1)
    
    ##Beta's part: uniform
    ll_b <- log(1)
    
    ##Phi's inverse-Wishart part: Phi
    ll_p <- log(diwish(matrix(phitrack[i,], ncol = d*q), v = v_2, S = Lambda_2))
    
    ##Put them all together and what do you get? The full posterior log likelihood
    full_post_log_like_track[i-1,] <- c(ll_wgy, ll_pgp, ll_m, ll_b, ll_p)
    
  }
  
  return(full_post_log_like_track)
  
}

#Trueloglikelihood (M-Logit)
multlogit_truloglike <- function(W, X, Z, randeff, y, mu, beta, psi, phi, test = FALSE){
  
  
  #Some dimensions
  p <- ncol(X) #number of covariates, including intercept
  d <- ncol(W) - 1 #number of multivariate categories
  N <- nrow(X) #number of total observations
  q <- randeff #number of random player effects (default is simply an intercept)
  m <- ncol(Z)/randeff #number of players
  
  #Algorithm uses a wide X
  X <- t(X) #should be p by N
  
  #Hyperparameters for prior distribution of Phi
  v_2 <- d*q
  Lambda_2 <- diag(d*q) + matrix(1, nrow = d*q)%*%matrix(1, ncol = d*q)
  
  #Multinomial part: W|Y
  ll_wgy <- 0
  for(l in 1:nrow(W)){
    sumexpyip1 <- 1 + sum(exp(y[l,]))
    y_prob <- unlist(c((exp(y[l,])/sumexpyip1), (1/sumexpyip1)))
    ll_wgyi <- dmultinom(W[l,], prob = y_prob, log = TRUE)
    ll_wgy <- ll_wgy + ll_wgyi
  }
  
  if(test == TRUE){
    browser()
  }
  
  ##Psi's normal part: Psi|Phi
  ll_pgp <- dmvn(as.vector(t(psi)), mu = rep(mu, m), sigma = kronecker(diag(m), matrix(as.vector(phi), byrow = F, ncol = d*q)), log = TRUE)
  
  ##Beta's part: Beta|Sigma (depends on block)
  #Block update Beta prior is uniform
  ll_m <- log(1)
  
  ll_b <- log(1)
  
  ##Phi's inverse-Wishart part: Phi
  ll_p <- log(diwish(matrix(as.vector(phi), ncol = d*q), v = v_2, S = Lambda_2))
  
  return(c(ll_wgy, ll_pgp, ll_m, ll_b, ll_p))
  
}

