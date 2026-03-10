
misd <- function(tpp, 
                 grange, 
                 tot_time,
                 num_iter, 
                 bw = 3, 
                 normN = T, 
                 verbose = T, 
                 makeplot = F) {
  supDist <- function (x, y) return (max (abs (x - y)))
  if (min(grange) > 0) grange <- c(0, grange)
  N <- length(tpp)
  delta_t <- diff(grange)  # bin width of histogram estimator g based on input grange
  g <- double(length(delta_t))  # intialize triggering function
  tdiffsmtx <- outer(tpp, tpp, FUN = "-")  # N by N matrix of times differences for tpp 
  tdiffs <- tdiffsmtx[lower.tri(tdiffsmtx, diag = F)]	
  
  A <- list(0)  
  bin_ind <- rep(length(A) + 1, length(tdiffs))
  converge = 0
  for(i in 2:length(grange)) {
    inds <- which(tdiffs > grange[i-1] & tdiffs <= grange[i]) # Find the indices in each bin
    if ((length(inds) == 1) & is.na(inds[1])) inds <- 0
    A[[i-1]] <- inds
    bin_ind[inds] <- i - 1
  }  
  
  # intitialize p_{ij} = 1/i for j <= i and 0 otherwise
  Pmtx_new <- matrix(0, nrow = N, ncol = N)
  for(i in 1:N) {
    Pmtx_new[i,] <- rep(1/i, N)
  }
  P_new <- Pmtx_new[lower.tri(Pmtx_new)] #note: diag = F is the default
  
  # start convergence algorithm 
  for(k in 1:num_iter) {
    Pmtx_old <- Pmtx_new
    P_old <- P_new		
    
    for(j in 1:length(delta_t)) {
      Z <- delta_t[j] * ifelse(normN, N, sum(P_old))
      g[j] <- sum(P_old[A[[j]]]) / Z   
    }
    
    
    mu <- sum(diag(Pmtx_old)) / tot_time
    
    # take the probabilities from the M-step and put them back in their right place in the g matrix using the bin indicator 
    fix <- g[bin_ind]
    fix[is.na(fix)] <- 0 
    
    gmtx <- matrix(0, nrow = N, ncol = N)
    gmtx[lower.tri(gmtx)] <- fix 
    
    # set the diagonal elements of the g matrix to the new background rate from the E-step 
    diag(gmtx) <- mu
    
    # E-step
    Pmtx_new <- matrix(0, nrow = N, ncol = N)  
    for(i in 1:N) {
      Pmtx_new[i,] <- gmtx[i,] / sum(gmtx[i,])
    }
    
    P_new <- Pmtx_new[lower.tri(Pmtx_new)]
    
    if(k > 2 & k %%10 == 0) {
      sup <- supDist(Pmtx_old, Pmtx_new)
      if(verbose) {
        cat (
          "Iteration: ", k,
          "SupDist: ", formatC(sup, digits = 8, width = 12, forma = "f"),
          "\n")
      }
      if( sup < 1e-5 ) {
        converge = 1
        break
      }	
    }
  }
  if (makeplot) {
    plot(grange[-1], g)
    lines(grange[-1], gauss.filt(g, bins = grange[-1], bw = bw)) 
  }
  return(list(mu = mu, g = g, grange = grange, delta_t = delta_t, converge = converge))	
  
}


#### NOT TESTED
misdtxy <- function(tpp, ppxy, granget, grangexy, distMat,
                    tot_time, num_iter, bw = 3, normN = T, verbose = T, plott = F) {
  supDist <- function (x, y) return (max (abs (x - y)))
  if (min(grange) > 0) granget <- c(0, granget)
  N <- length(tpp)
  delta_t <- diff(granget)  # bin width of histogram estimator g based on input grange
  delta_t <- diff(grangexy)
  
  g <- double(length(delta_t))  # intialize triggering function
  tdiffsmtx <- outer(tpp, tpp, FUN = "-")  # N by N matrix of times differences for tpp 
  tdiffs <- tdiffsmtx[lower.tri(tdiffsmtx, diag = F)]	
  
  xydiffs <- distMat[lower.tri(distMat)]
  
  A <- list(0)  # for each bin, stores vector of event ids
  bin_ind <- rep(length(A) + 1, length(tdiffs)) # for each event, labels which bin it belongs in
  converge = 0
  for(i in 2:length(grange)) {
    inds <- which(tdiffs > grange[i-1] & tdiffs <= grange[i]) # Find the indices in each bin
    if ((length(inds) == 1) & is.na(inds[1])) inds <- 0
    A[[i-1]] <- inds
    bin_ind[inds] <- i - 1
  }  
  
  # intitialize p_{ij} = 1/i for j <= i and 0 otherwise
  Pmtx_new <- matrix(0, nrow = N, ncol = N)
  for(i in 1:N) {
    Pmtx_new[i,] <- rep(1/i, N)
  }
  P_new <- Pmtx_new[lower.tri(Pmtx_new)] #note: diag = F is the default
  
  # start convergence algorithm 
  for(k in 1:num_iter) {
    Pmtx_old <- Pmtx_new
    P_old <- P_new		
    
    for(j in 1:length(delta_t)) {
      Z <- delta_t[j] * ifelse(normN, N, sum(P_old))
      g[j] <- sum(P_old[A[[j]]]) / Z   
    }
    
    
    mu <- sum(diag(Pmtx_old)) / tot_time
    
    # take the probabilities from the M-step and put them back in their right place in the g matrix using the bin indicator 
    fix <- g[bin_ind]
    fix[is.na(fix)] <- 0 
    
    gmtx <- matrix(0, nrow = N, ncol = N)
    gmtx[lower.tri(gmtx)] <- fix 
    
    # set the diagonal elements of the g matrix to the new background rate from the E-step 
    diag(gmtx) <- mu
    
    # E-step
    Pmtx_new <- matrix(0, nrow = N, ncol = N)  
    for(i in 1:N) {
      Pmtx_new[i,] <- gmtx[i,] / sum(gmtx[i,])
    }
    
    P_new <- Pmtx_new[lower.tri(Pmtx_new)]
    
    if(k > 2 & k %%10 == 0) {
      sup <- supDist(Pmtx_old, Pmtx_new)
      if(verbose) {
        cat (
          "Iteration: ", k,
          "SupDist: ", formatC(sup, digits = 8, width = 12, forma = "f"),
          "\n")
      }
      if( sup < 1e-5 ) {
        converge = 1
        break
      }	
    }
  }
  if (plott) {
    plot(grange[-1], g)
    lines(grange[-1], gauss.filt(g, bins = grange[-1], bw = bw)) 
  }
  return(list(mu = mu, g = g, grange = grange, delta_t = delta_t, converge = converge))	
  
}

