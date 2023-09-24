#Prediction functions for mmln output

newobsrandeff <- function(np_W, np_X, np_Z, randeff, base, iter, mhiter, mhscale, proposal = c("normal", "beta", "normbeta"), block = TRUE, sigma_chain, phi_chain, beta_chain, trY = NULL, trPsi = NULL, fixY = FALSE, fixPsi = FALSE, whichYstart = c("mean", "inity"), whichCOVprior = c("uninf", "weak"), test = FALSE, oldway = FALSE){
  
  if(oldway == TRUE){
    
    np_Psi_mat <- vector("list", nrow(sigma_chain))
    for(i in 1:nrow(sigma_chain)){
      temp_fixSigma <- matrix(sigma_chain[i,], ncol = ncol(np_W) - 1, byrow = T)
      temp_fixPhi <- matrix(phi_chain[i,], ncol = ncol(np_W) - 1, byrow = T)
      temp_fixBeta <- matrix(beta_chain[i,], ncol = ncol(np_W) - 1, byrow = T)
      
      temp_res.estphi <- mixgibbs(W = np_W, X = np_X, Z = np_Z, randeff = randeff, base = base, iter = iter, mhiter = mhiter, mhscale = mhscale, proposal = proposal, block = block, r.seed = i, trY = trY, trBeta = t(temp_fixBeta), trPsi = t(t(trPsi)), trSigma = temp_fixSigma, trPhi = temp_fixPhi, startru = F, fixY = fixY, fixBeta = T, fixPsi = fixPsi, fixSigma = T, fixPhi = T, whichYstart = whichYstart, test = test)
      np_Psi_mat[[i]] <- temp_res.estphi$Psi[iter + 1,]
      
    }
    
    np_Psi_mat <- do.call(rbind, np_Psi_mat)
    return(np_Psi_mat)
    
  } else{
    
    temp_fixSigma <- matrix(colMeans(sigma_chain), ncol = ncol(np_W) - 1, byrow = T)
    temp_fixBeta <- matrix(colMeans(beta_chain), ncol = ncol(np_W) - 1, byrow = T)
    temp_fixPhi <- matrix(colMeans(phi_chain), ncol = ncol(np_W) - 1, byrow = T)
    
    temp_res.estphi <- mixgibbs(W = np_W, X = np_X, Z = np_Z, randeff = randeff, base = base, iter = iter, mhiter = mhiter, mhscale = mhscale, proposal = proposal, block = block, trY = trY, trBeta = t(temp_fixBeta), trPsi = trPsi, trPhi = temp_fixPhi, trSigma = temp_fixSigma, startru = F, fixY = fixY, fixSigma = T, fixBeta = T, fixPsi = fixPsi, fixPhi = T, whichYstart = whichYstart, whichCOVprior = whichCOVprior, test = test)
    np_Psi_mat <- temp_res.estphi$Psi
    
    return(np_Psi_mat)
    
  }
}


#Mult-Logit Random Effect OOS Recovery
newobsrandeff.multlogit <- function(np_W, np_X, np_Z, randeff, base, iter, mhiter_psi, mhscale_psi, phi_chain, beta_chain, mu_chain, trY = NULL, trPsi = NULL, fixY = FALSE, fixPsi = FALSE, whichYstart = c("mean", "inity"), test = FALSE, oldway = FALSE){
  
  if(test == TRUE){ browser() }
  
  if(oldway == TRUE){
    
    np_Psi_mat <- vector("list", nrow(phi_chain))
    for(i in 1:nrow(phi_chain)){
      temp_fixPhi <- matrix(phi_chain[i,], ncol = ncol(np_W) - 1, byrow = T)
      temp_fixBeta <- rbind(matrix(mu_chain[i,], ncol = ncol(np_W) - 1, byrow = T), matrix(beta_chain[i,], ncol = ncol(np_W) - 1, byrow = T))
      
      temp_res.estphi <- multmixgibbs_b0prior(W = np_W, X = np_X, Z = np_Z, randeff = randeff, base = base, iter = iter, mhiter_psi = mhiter_psi, mhscale_psi = mhscale_psi, trY = trY, trBeta = t(temp_fixBeta), trPsi = trPsi, trPhi = temp_fixPhi, fixIntercept = T, fixBeta = T, fixPsi = fixPsi, fixPhi = T, whichYstart = whichYstart, test = test)
      np_Psi_mat[[i]] <- temp_res.estphi$Psi[iter + 1,]
      
    }
    
    np_Psi_mat <- do.call(rbind, np_Psi_mat)
    return(np_Psi_mat)
    
  } else{
    
    temp_fixBeta <- rbind(colMeans(mu_chain), matrix(colMeans(beta_chain), ncol = ncol(np_W) - 1, byrow = T))
    temp_fixPhi <- matrix(colMeans(phi_chain), ncol = ncol(np_W) - 1, byrow = T)
    
    temp_res.estphi <- multmixgibbs_b0prior(W = np_W, X = np_X, Z = np_Z, randeff = randeff, base = base, iter = iter, mhiter_psi = mhiter_psi, mhscale_psi = mhscale_psi, trY = trY, trBeta = t(temp_fixBeta), trPsi = trPsi, trPhi = temp_fixPhi, fixIntercept = T, fixBeta = T, fixPsi = fixPsi, fixPhi = T, whichYstart = whichYstart, test = test)
    np_Psi_mat <- temp_res.estphi$Psi
    
    return(np_Psi_mat)
    
  }
}


#IN THE BELOW, THE np_pi_pred will be wrong if there's a base other than the last category
predict.randeffseparate <- function(np_pa, np_X, np_Z, dimW, sigma_chain, phi_chain, beta_chain, psi_chain, cycle = FALSE, mean = FALSE, marginal = FALSE, singleiter = FALSE, whichiter = 500, test = FALSE){
  
  if(test == TRUE){browser()}
  
  y_ij <- NULL
  
  no.samples <- nrow(sigma_chain)
  
  if(singleiter == TRUE){
    y_ij <- rmvn(n = no.samples, mu = t(np_X)%*%matrix(beta_chain[whichiter,], ncol = dimW - 1, byrow = T) + psi_chain[whichiter,], sigma = matrix(sigma_chain[whichiter,], ncol = dimW - 1, byrow = T))
    pi_ij <- t(apply(y_ij, 1, ytopi))
    w_ij <- t(apply(pi_ij, 1, rmultinom, n = 1, size = np_pa))
    if(marginal == TRUE){
      sigma <- matrix(sigma_chain[whichiter,], ncol = dimW - 1, byrow = T)
      phi <- matrix(phi_chain[whichiter,], ncol = dimW - 1, byrow = T)
      y_ij <- rmvn(n = no.samples, mu = t(np_X)%*%matrix(beta_chain[whichiter,], ncol = dimW - 1, byrow = T), sigma = sigma + phi)
      pi_ij <- t(apply(y_ij, 1, ytopi))
      w_ij <- t(apply(pi_ij, 1, rmultinom, n = 1, size = np_pa))
    }
  }
  
  if(mean == TRUE){
    postmean.sigma <- matrix(colMeans(sigma_chain), ncol = dimW - 1, byrow = T)
    postmean.beta <- matrix(colMeans(beta_chain), ncol = dimW - 1, byrow = T)
    postmean.psi <- matrix(colMeans(psi_chain), ncol = dimW - 1)
    y_ij <- rmvn(n = no.samples, mu = t(np_X)%*%postmean.beta + postmean.psi, sigma = postmean.sigma)
    pi_ij <- t(apply(y_ij, 1, ytopi))
    w_ij <- t(apply(pi_ij, 1, rmultinom, n = 1, size = np_pa))
    if(marginal == TRUE){
      postmean.phi <- matrix(colMeans(phi_chain), ncol = dimW-1, byrow = T)
      y_ij <- rmvn(n = no.samples, mu = t(np_X)%*%postmean.beta, sigma = postmean.phi)
      pi_ij <- t(apply(y_ij, 1, ytopi))
      w_ij <- t(apply(pi_ij, 1, rmultinom, n = 1, size = np_pa))
    }
  }
  
  if(cycle == TRUE){
    mu_ij <- t(apply(beta_chain, 1, function(x){ t(np_X)%*%matrix(x, ncol = dimW - 1, byrow = T) })) + psi_chain
    epsilon_ij <- t(apply(sigma_chain, 1, function(x){ rmvn(n = 1, mu = rep(0, dimW - 1), sigma = matrix(x, ncol = dimW - 1, byrow = T)) }))
    y_ij <- mu_ij + epsilon_ij
    pi_ij <- t(apply(y_ij, 1, ytopi))
    w_ij <- t(apply(pi_ij, 1, rmultinom, n = 1, size = np_pa))
    if(marginal == TRUE){
      psi_chain <- t(apply(phi_chain, 1, function(x){ rmvn(n = 1, mu = rep(0, dimW-1), sigma = matrix(x, ncol = dimW - 1, byrow = T))}))
      y_ij <- t(apply(beta_chain, 1, function(x){ t(np_X)%*%matrix(x, ncol = dimW - 1, byrow = T) })) + psi_chain + epsilon_ij
      pi_ij <- t(apply(y_ij, 1, ytopi))
      w_ij <- t(apply(pi_ij, 1, rmultinom, n = 1, size = np_pa))
    }
  }
  
  if(is.null(y_ij)){return("No Scheme Specified")} else{
    result <- list(predy = y_ij, predpi = pi_ij, predw = w_ij, psisamples = psi_chain)
    return(result)
  }
  
}

#Mult-Logit Prediction
predict.randeffseparate.multlogit <- function(np_pa, np_X, np_Z, dimW, phi_chain, mu_chain, beta_chain, psi_chain, cycle = FALSE, mean = FALSE, marginal = FALSE, singleiter = FALSE, whichiter = 500, test = FALSE, no.samples = 1000){
  
  if(test == TRUE){browser()}
  
  y_ij <- NULL
  
  if(singleiter == TRUE){
    y_ij <- t(np_X)[,-1]%*%matrix(beta_chain[whichiter,], ncol = dimW - 1, byrow = T) + psi_chain[whichiter,]
    pi_ij <- ytopi(y_ij)
    w_ij <- t(rmultinom(n = no.samples, size = np_pa, prob = pi_ij))
    if(marginal == TRUE){
      phi <- matrix(phi_chain[whichiter,], ncol = dimW-1, byrow = T)
      y_ij <- rmvn(n = no.samples, mu = t(np_X)[,-1]%*%matrix(beta_chain[whichiter,] + mu_chain[whichiter,], ncol = dimW - 1, byrow = T), sigma = phi)
      pi_ij <- t(apply(y_ij, 1, ytopi))
      w_ij <- t(apply(pi_ij, 1, rmultinom, n = 1, size = np_pa))
    }
  }
  
  if(mean == TRUE){
    postmean.beta <- matrix(colMeans(beta_chain), ncol = dimW - 1, byrow = T)
    postmean.psi <- matrix(colMeans(psi_chain), ncol = dimW - 1)
    y_ij <- t(np_X)[,-1]%*%postmean.beta + postmean.psi
    pi_ij <- ytopi(y_ij)
    w_ij <- t(rmultinom(n = no.samples, size = np_pa, prob = pi_ij))
    if(marginal == TRUE){
      postmean.phi <- matrix(colMeans(phi_chain), ncol = dimW-1, byrow = T)
      y_ij <- rmvn(n = no.samples, mu = t(np_X)[,-1]%*%postmean.beta + colMeans(mu_chain), sigma = postmean.phi)
      pi_ij <- t(apply(y_ij, 1, ytopi))
      w_ij <- t(apply(pi_ij, 1, rmultinom, n = 1, size = np_pa))
    }
  }
  
  if(cycle == TRUE){
    y_ij <- t(apply(beta_chain, 1, function(x){ t(np_X)[,-1]%*%matrix(x, ncol = dimW - 1, byrow = T) })) + psi_chain
    pi_ij <- t(apply(y_ij, 1, ytopi))
    w_ij <- t(apply(pi_ij, 1, rmultinom, n = 1, size = np_pa))
    if(marginal == TRUE){
      beta_chain <- cbind(mu_chain, beta_chain)
      psi_chain <- t(apply(phi_chain, 1, function(x){ rmvn(n = 1, mu = rep(0, dimW - 1), sigma = matrix(x, ncol = dimW - 1, byrow = T))}))
      y_ij <- t(apply(beta_chain, 1, function(x){ t(np_X)%*%matrix(x, ncol = dimW - 1, byrow = T) })) + psi_chain
      pi_ij <- t(apply(y_ij, 1, ytopi))
      w_ij <- t(apply(pi_ij, 1, rmultinom, n = 1, size = np_pa))
      psi_chain <- mu_chain + psi_chain
    }
  }
  
  if(is.null(y_ij)){return("No Scheme Specified")} else{
    result <- list(predy = y_ij, predpi = pi_ij, predw = w_ij, psisamples = psi_chain)
    return(result)
  }
  
}