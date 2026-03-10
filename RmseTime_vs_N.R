# Compare time/RMSE of unregularized analytic, empirical and my regularized estimator
# for variable productivity with fixed br ratio

source("simetas.R")
source("Kestimate_functions.R")
library(RSpectra)

######### Generate 
set.seed(375) 
# ba <- 140; bc = 0.7 * ba
theta_beta = 5
a3 = 5; b3 = 2; # params for gaussian triggering
mu = 1; theta_K = 0.5; 
BigT = 1000

test.times <- seq(0.8, 2.2, by = 0.05) * 1000
FixedBr <- 0.1
MyNs <- c()
Brs <- c()
todrop <- c()
nreps <- 2 # for each size, generate nreps simulations

SvdRmse_avg <- c() # Average rmse over nreps
SolveRmse_avg <- c()
EmpRmse_avg <- c()
OlsRmse_avg <- c() # small delta t
OlsRmse_avg2 <- c() # big(ger) delta t
UnregRmse_avg <- c()
UnregRmse_avg2 <- c() # unscaled unregularized


SvdTimes_avg <- c()
SolveTimes_avg <- c()
EmpTimes_avg <- c()
OlsTimes_avg <- c() # small delta t
OlsTimes_avg2 <- c() # big(ger) delta t
UnregTimes_avg <- c()
EmpTimes_avg <- c()

MeanNs <- c() # mean for each BigT
MeanBrs <- c()
UnregFails <- double(length(test.times))

for (tt in 1:length(test.times)) {
  BigT <- test.times[tt]
  print(paste(tt, "/", length(test.times)))
  
  xs <- seq(0, BigT, length = 10000) # interval to interpolate estimates to for comparability
  SvdAvg <- double(length(xs))
  EmpAvg <- double(length(xs))
  SolveAvg <- double(length(xs))
  OlsAvg <- double(length(xs))
  UnregAvg <- double(length(xs))
  
  SvdRmse <- c()
  SolveRmse <- c()
  EmpRmse <- c()
  OlsRmse <- c()
  UnregRmse <- c()
  UnregRmse2 <- c() # unscaled unregularized
  
  SvdTimes <- c()
  SolveTimes <- c()
  EmpTimes <- c()
  OlsTimes <- c()
  UnregTimes <- c()
  EmpTimes <- c()
  
  zns <- c()
  ## K \sim ba * N(bb, be) + bc * N(bd, bf)
  ## need to adjust as you increase magnitude so it doesn't become explosive
  bb = 0.25 * BigT; bd = 0.75 * BigT; be = 0.07 * BigT; bf = be;
  ba <- FixedBr / mean(dnorm(xs, m = bb, s = be) + 0.5 * dnorm(xs, m = bd, s = bf))
  bc = 0.5 * ba
  myktest = ba * dnorm(xs, mean=bb, sd=be) + bc * dnorm(xs, mean=bd, sd=bf)
  
  if (max(myktest) > 1) {
    print("Done")
    break
  }
  BW <- BigT / 10
  MeanBrs <- c(MeanBrs, mean(myktest))
  
  for (rep in 1:nreps) {
  z0 = simhawk(BigT = BigT, gmi = normprod, gxy = pointxy, gt = expgt, 
               theta = list(mu = mu, K = theta_K, beta = theta_beta, b = 5,
                            ba=ba, bb=bb, bc=bc, bd=bd, be=be, bf=bf) )   
  n <- z0$n
  zns <- c(zns, n)
  t <- z0$t
  myk = ba * dnorm(t, mean=bb, sd=be) + bc * dnorm(t, mean=bd, sd=bf)
  
  ##############################
  #### Pointwise estimates with SVD
  start <- Sys.time()
  G = matrix(0, nrow = n - 1, ncol = n - 1) ##G[1,1] = g(t2-t1). deriv wrt beta1 for 1/lam2, 1/lam3.
  for (i in 1:(n-1)) {
    for (j in i:(n - 1)) {
      G[i,j] = dexp(t[j + 1] - t[i], rate = 1/theta_beta)
    }}
  Gs <- svd(G)

  # Select optimal regularization param (increase until sum > n - muT)
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
  SvdTimes <- c(SvdTimes, difftime(Sys.time(), start, units = "secs"))
  SvdAvg <- SvdAvg + ksmooth(t, Kest, kernel = "normal", b = BW, x.points = xs)$y
  SvdRmse <- c(SvdRmse, mse(myk - Kest))

  #######################
  #### Empirical estimator
  start <- Sys.time()
  incr <- 7
  Kemp <- rep(0, n)
  for(i in 1:n){ 
    Kemp[i] = max(0, sum((t > t[i]) & (t < t[i] + incr)) - incr * mu) ## how many extra points occurred
  } 
  Kemp <- (Kemp / sum(Kemp)) * (n - mu * BigT)
  EmpAvg <- EmpAvg + ksmooth(t, Kemp, kernel = "normal", b = BigT / 10, x.points = xs)$y
  EmpTimes <- c(EmpTimes, difftime(Sys.time(), start, units = "secs"))
  Kemp <- ksmooth(t, Kemp, kernel = "normal", b = BW, x.points = t)$y
  EmpRmse <- c(EmpRmse, mse(myk - Kemp))
  
  ###################################
  ## Unregularized solve
  start <- Sys.time()
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
    if (sum(Kest4) > 1e-1) {
      UnregRmse2 <- c(UnregRmse2, mse(myk - Kest4))
      Kest4 <- (Kest4 / sum(Kest4)) * (n - mu * BigT)
      UnregAvg <- UnregAvg + ksmooth(t, Kest4, kernel = "normal", b = BW, x.points = xs)$y
      UnregTimes <- c(UnregTimes, difftime(Sys.time(), start, units = "secs"))
      UnregRmse <- c(UnregRmse, mse(myk - Kest4))
    }
  }
  
  
  ###############################
  ############### Binned OLS: small discretization
  start <- Sys.time()
  I <- 10
  delt <- BigT * (1.5 / 1000); # discretization parameters
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
  OlsAvg <- OlsAvg + ksmooth(t, Kest3, b = BW, x.points = xs)$y
  OlsTimes <- c(OlsTimes, difftime(Sys.time(), start, units = "secs"))
  OlsRmse <- c(OlsRmse, mse(myk - Kest3))
  }    
  
  MeanNs <- c(MeanNs, mean(zns))
  
  SvdTimes_avg <- c(SvdTimes_avg, mean(SvdTimes))
  EmpTimes_avg <- c(EmpTimes_avg, mean(EmpTimes))
  OlsTimes_avg <- c(OlsTimes_avg, mean(OlsTimes))
  UnregTimes_avg <- c(UnregTimes_avg, mean(UnregTimes))
  
  SvdRmse_avg <- c(SvdRmse_avg, mean(SvdRmse))
  EmpRmse_avg <- c(EmpRmse_avg, mean(EmpRmse))
  OlsRmse_avg <- c(OlsRmse_avg, mean(OlsRmse))
  UnregRmse_avg <- c(UnregRmse_avg, mean(UnregRmse))
  UnregRmse_avg2 <- c(UnregRmse_avg2, mean(UnregRmse2))
  if (length(UnregRmse) < nreps) UnregFails[tt] <- mean(zns)
}



######
par(mfrow = c(1, 2))
## Time vs n compare
plot(MeanNs, SvdTimes_avg, ylim = c(0, max(SvdTimes_avg, OlsTimes_avg)),
     main = paste("Mean Productivity", FixedBr), cex.axis = 1.3, cex.lab=1.2,
     xlab = "n", ylab = "Estimation Time (s)", pch = 16, lwd=2, col = "#1f77b4")
lines(MeanNs, SvdTimes_avg, col = "#1f77b4")

# points(MeanNs, UnregTimes_avg, pch = 2)
points(MeanNs, UnregTimes_avg, pch = 0, lwd=2, col = "#ff7f0e")
lines(MeanNs, UnregTimes_avg, pch = 0, col = "#ff7f0e")

points(MeanNs, EmpTimes_avg, pch = 2,  lwd=2, col = "#2ca02c")
lines(MeanNs, EmpTimes_avg, pch = 2,  col = "#2ca02c")

# points(MeanNs, SolveTimes_avg, pch = 5)
points(MeanNs + 10, OlsTimes_avg, pch = 5,  lwd=2, col = "#d62728")
lines(MeanNs + 10, OlsTimes_avg, pch = 5, col = "#d62728")

# points(MeanNs + 10, OlsTimes_avg2, pch = 5, col = "blue")
legend("topleft", pch = c(16, 0, 2, 5), cex = 1.2, lwd=2, lty=0, bty="n",
       col = c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728"),
       legend = c("Regularized", "Unregularized",
                  "Empirical", "Binned Least-Squares"))

## Rmse vs n compare
plot(MeanNs, SvdRmse_avg, pch = 16, lwd=2, col = "#1f77b4",
     ylim = range(c(UnregRmse_avg, SvdRmse_avg, EmpRmse_avg, OlsRmse_avg)),
     xlab = "n", ylab = "RMSE", cex.axis = 1.2, cex.lab = 1.2)
lines(MeanNs, SvdRmse_avg, col = "#1f77b4")

points(MeanNs, UnregRmse_avg, pch = 0, lwd=2, col = "#ff7f0e")
lines(MeanNs, UnregRmse_avg, pch = 0, col = "#ff7f0e")

points(MeanNs, EmpRmse_avg, pch = 2, lwd=2, col = "#2ca02c")
lines(MeanNs, EmpRmse_avg, pch = 2, col = "#2ca02c")

points(MeanNs, OlsRmse_avg, pch = 5, lwd=2, col = "#d62728")
lines(MeanNs, OlsRmse_avg, pch = 5, col = "#d62728")
