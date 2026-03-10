### Generate 100 simulations, fit empirical/unreg/regularized/least squares estimator to
### each and plot average estimates.

source("Kestimate_functions.R")
source("simetas.R")
library(RSpectra)

######### Generate 
set.seed(375) 
theta_beta = 5
a3 = 5; b3 = 2; # params for gaussian triggering
mu = 1; theta_K = 0.25; 
bm = 300; bs = 200;
BigT = 1000

####### Repeat estimation M times
M <- 50
zns <- c()

xs <- seq(0, BigT, length = 10000) # interval to interpolate estimates to for comparability
SvdAvg <- double(length(xs))
EmpAvg <- double(length(xs))
SolveAvg <- double(length(xs))
OlsAvg <- double(length(xs))
OlsAvg2 <- double(length(xs))
UnregAvg <- double(length(xs))

SvdTimes <- c()
SolveTimes <- c()
OlsTimes <- c()
UnregTimes <- c()

SvdRmse <- c()
EmpRmse <- c()
SolveRmse <- c()
OlsRmse <- c()
UnregRmse <- c()

singulars <- 0 # for solve
rparams <- c() # for svd

for (rep in 1:M) {
  if ((rep %% 10) == 0) print(rep)
  
  ## normprod: K \sim ba * N(bb, be) + bc * N(bd, bf)
  bb = 0.25*BigT; bd = 0.75*BigT; 
  be = 0.07 * BigT; bf = be;
  ba <- 0.7 / max(dnorm(xs, m = bb, s = be) + 0.5 * dnorm(xs, m = bd, s = bf))
  bc = 0.5 * ba
  
  ## normprod2: K \sim bw * N(bm, bs) 
  bw = 0.5 / max(dnorm(xs, m = bm, s = bs))
  
  
  z0 = simhawk(BigT = BigT, gmi = normprod2, gxy = pointxy, gt = expgt, 
               theta = list(mu = mu, K = theta_K, beta = theta_beta, b = 5,
                            ba=ba, bb=bb, bc=bc, bd=bd, be=be, bf=bf,
                            bw=bw, bm=bm, bs=bs) )   
  n <- z0$n
  zns <- c(zns, n)
  BW = BigT / 10
  
  ##################
  #### SVD
  start <- Sys.time()
  t <- z0$t
  G = matrix(0, nrow = n - 1, ncol = n - 1) ##G[1,1] = g(t2-t1). deriv wrt beta1 for 1/lam2, 1/lam3. 
  for (i in 1:(n-1)) {
    for (j in i:(n - 1)) {
      G[i,j] = dexp(t[j+1] - t[i], rate = 1/theta_beta)
    }} 
  Gs <- svd((G))
  
  # Select optimal regularization param (achieves sufficiently large K)
  rtests <- seq(0, 0.3, by = 0.01)[-1]
  reg_params <- rtests * Gs$d[1]

  Kprev <- c()
  for (ind in 1:length(reg_params)) {
    reg <- reg_params[ind]^2
    
    inlam0 <- tikh(Gs$u, Gs$d, Gs$v, rep(1, n - 1), reg)
    lambda0 <- c(mu, drop(1/inlam0$coef))
    Kest0 <- tikh(Gs$v, Gs$d, Gs$u, lambda0[2:n] - lambda0[1], reg)
    Kest0 <- c(drop(Kest0$coef), 0)

    Kest <- ksmooth(t, Kest0, kernel = "normal", b = BW, x.points = t)$y
    Kest <- force0(Kest)
    if (sum(Kest) > (n - mu * BigT)) break
    else Kprev <- Kest
  }
  if (ind > 1) Kest <- Kprev
  rparams <- c(rparams, rtests[ind - 1])
  SvdTimes <- c(SvdTimes, difftime(Sys.time(), start, units = "secs"))
  SvdAvg <- SvdAvg + ksmooth(t, Kest, kernel = "normal", b = BW, x.points = xs)$y
  
  ################################
  ## Empirical estimator
  incr <- 7
  Kemp <- rep(0, n)
  for(i in 1:n){ 
    Kemp[i] = max(0, sum((t > t[i]) & (t < t[i] + incr)) - incr * mu) ## how many extra points occurred
  } 
  Kemp <- (Kemp / sum(Kemp)) * (n - mu * BigT)
  EmpAvg <- EmpAvg + ksmooth(t, Kemp, kernel = "normal", b = BW, x.points = xs)$y
  
  ####################################
  ## Unregularized solve
  start <- Sys.time()
  t <- z0$t
  G = matrix(0, nrow = n - 1, ncol = n - 1) ##G[1,1] = g(t2-t1). deriv wrt beta1 for 1/lam2, 1/lam3.
  for (i in 1:(n-1)) {
    for (j in i:(n - 1)) {
      G[i,j] = dexp(t[j+1] - t[i], rate = 1/theta_beta)
    }}
  if (rcond(G) == 0) {
    warning(paste("Singular matrix with BigT=", BigT))
  } else {
    inlam2 <- solve(G, rep(1, n - 1))
    lambda2 <- c(mu, 1/inlam2)
    Kest4 <- solve(t(G),  lambda2[2:n] - lambda2[1])
    Kest4 <- ksmooth(t, c(Kest4, 0), kernel = "normal", b = BW, x.points = t)$y
    Kest4 <- force0(Kest4)
    if (sum(Kest4) > 0) {
      Kest4 <- (Kest4 / sum(Kest4)) * (n - mu * BigT)
      UnregAvg <- UnregAvg + ksmooth(t, Kest4, kernel = "normal", b = BW, x.points = xs)$y
      UnregTimes <- c(UnregTimes, difftime(Sys.time(), start, units = "secs")) 
    } else {
      singulars <- singulars + 1
    }
  }
  
  ###############################
  ############### Binned OLS: coarse discretization
  I <- 10; delt <- 1.5; # discretization parameters
  bin.breaks <- seq(0, BigT + delt, by = delt)
  T2 <- length(bin.breaks) - 1
  W <- ceiling(T2 / I)
  true.ts = bin.breaks[seq(1, T2, by = I)] # for smoothing; the times each interval corresponds to
  
  hout <- hist(z0$t, breaks = bin.breaks, plot = F)
  binCounts <- hout$counts
  
  G <- matrix(0, nrow = T2 - 1, ncol = W)
  for (i in 2:T2) {
    upTo <- min(W, ceiling(i / I)) # which period the day is in
    for(j in 1:upTo) { # which periods precede bin i
      daysInJ <- seq((j-1) * I + 1, min(i - 1, j * I))
      daysInJ <- daysInJ[daysInJ < i] # some edge cases allow this
      G[i - 1, j] <- sum(dexp(i - daysInJ, 1/theta_beta) * binCounts[daysInJ])
    }
  }
  
  sG <- svds(G, k = 1, nu = 0, nv = 0)$d # max singular value
  reg_param <- 0.01 * sG
  L <- reg_param * diag(1, ncol = ncol(G), nrow = ncol(G))
  
  Kest0 <- (qr.solve(rbind(G, L), c(binCounts[-1] - mu * delt, rep(0, nrow(L)))))
  Kest0 <- c(Kest0, 0)
  Kest3 <- ksmooth(true.ts, Kest0, kernel = "normal", b = BW, x.points = t)$y
  Kest3 <- force0(Kest3)
  OlsAvg <- OlsAvg + ksmooth(t, Kest3, kernel = "normal", b = BW, x.points = xs)$y

}    

## plot
# myk <- ba * dnorm(xs, mean=bb, sd=be) + bc * dnorm(xs, mean=bd, sd=bf)
myk <- bw * dnorm(xs, m = bm, s = bs)
# myk <-rep(theta_K, length(xs))
mean(myk)

par(mfrow=c(1,1))
plot(xs, myk, type = "l", xlab = "t", ylab = "K(t)", lwd = 3,
     ylim = c(0, max(SvdAvg / M, EmpAvg / M, UnregAvg / M, myk)),
     # main = paste0("K = ", theta_K ))
     # main = paste0("K(t) ~ N(",bb,",",be,") + N(",bd,",",bf, ")"))
     main = paste0("K(t) ~ N(", bm, ",", bs, ")"))
lines(xs, SvdAvg / M, lty = 1, lwd=2, col = "#1f77b4")
lines(xs, UnregAvg / M, lty = 1, lwd=2, col = "#ff7f0e")
lines(xs + delt*I, OlsAvg / M, lty = 1, lwd=2, col = "#d62728")
lines(xs, EmpAvg / M, lty = 1, lwd=2, col = "#2ca02c")
legend("topright", lty = c(1, 1,1,1,1), 
       col = c(1, "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728"),
       bty="n", lwd=2, cex=1.2,
       legend = c("True K", "Regularized", "Unregularized", "Empirical", "Binned Least-Squares"))

