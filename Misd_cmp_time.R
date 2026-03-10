

source("simetas.r")
source("misd.R")
source("Kestimate_functions.R")
library(RSpectra)

#### generate 
X1 = 1; Y1 = 1
a3 = 7; b3 = 2 # gaussian triggering parameters
theta_beta = 5; # exp triggering parameters

## K \sim ba * N(bb, be) + bc * N(bd, bf)
ba = 80; bc = 40; # scaling
bb = 90; bd = 270 # mean
be = 21; bf = 21  # sd

mu = 1; theta_K = 0.5;
BigT = 1000

##########################################
#### Generate datasets at different horizons, increasing n
##########################################
test.times <- seq(0.8, 2, by = 0.1) * 1000
nt <- length(test.times)

zns <- double(nt)
grmses <- double(nt)
misd.grmses <- double(nt)
mu.err <- double(nt)
misd.mu.err <- double(nt)


ztimes <- double(nt)
misd.times <- double(nt)

for (rep in 1:nt) {
  
  BigT <- test.times[rep]
  print(BigT)
  
  # Binning parameters:
  t.delta <- 0.5
  day.breaks <- seq(0, BigT + t.delta, by = t.delta)
  ndays <- length(day.breaks) - 1
  W2.true <- 20 # the cont. interval you want to measure g on
  W2 <- W2.true / t.delta # on many discrete time intervals this represents
  true.ts <- day.breaks[1:W2 + 1] # the continuous times each interval represents
 
  misd.nbins <- 15
  grange <- seq(0, W2.true, length = misd.nbins + 1)
  delt <- diff(grange)[1] 
  
  # Generate
  z0 = simhawk(BigT = BigT, gmi = pointprod, gxy = pointxy, gt = normg, 
               theta = list(mu = mu, 
                            b = 5, # for mdensity
                            K = theta_K, # for uniform productivity
                            beta = theta_beta, # exp triggering
                            a3 = a3, b3 = b3 # gaussian triggering
               ) )
  n <- z0$n
  zns[rep] <- n
  bout <- hist(z0$t, day.breaks, plot = FALSE) # bin unique event times into "days"
  dailyCounts <- bout$counts
  
  #### Get Least-squares estimates
  start <- Sys.time()
  A <- matrix(0, nrow = ndays - 1, ncol = W2)
  for (i in 2:ndays) { 
    upTo <- min(i - 1, W2) # how many days before
    for(j in 1:upTo) { 
      dayInd <- i - j
      weekInd <- ceiling(dayInd / I) # which period the day is in
      A[i - 1, j] <- theta_K * dailyCounts[dayInd]
    }
  }
  A <- cbind(1, (A))
  
  sA <- 0.1 * svds(A, k = 1, nu = 0, nv = 0)$d # max sing val of A
  rparams[rep] <- sA
  
  L <- diag(c(0, rep(1, W2)), ncol = W2 + 1) # ridge
  gest0 <- qr.solve(rbind(A, sA * L), c(dailyCounts[-1], rep(0, nrow(L))) )
  
  mu.err[rep] <- norm(gest0[1] - mu)
  
  gest <- ksmooth(true.ts, gest0[-1], b = bw.nrd(true.ts), kernel = "normal", x.points = true.ts)$y
  if (all(gest <= 0)) stop("all negative")
  gest <- (gest) / sum(t.delta * force0(gest))
  
  gtrue <- dnorm(true.ts, m = a3, s = b3) / sum(t.delta * dnorm(true.ts, m = a3, s = b3))
  # gtrue <- dexp(true.ts, 1/theta_beta) / sum(t.delta * dexp(true.ts, 1/theta_beta))
  grmses[rep] <- norm(gtrue - gest)
  
  ztimes[rep] <- difftime(Sys.time(), start, units = "secs")
  
  #### Compare to MISD
  misd.start <- Sys.time()
  misd.out <- misd(z0$t, grange = grange, num_iter = 1e3, tot_time = BigT, verbose = F)
  misd.gest <- misd.out$g / sum(delt * misd.out$g)
  
  misd.gtrue <- dnorm(misd.out$g, m = a3, s = b3) / sum(delt * dnorm(misd.out$g, m = a3, s = b3))
  misd.grmses[rep] <- norm(misd.gtrue - misd.gest)
  misd.mu.err[rep] <- norm(misd.out$mu - mu)
  misd.times[rep] <- difftime(Sys.time(), misd.start, units = "secs")
}

res <- data.frame(list(n = zns, lsq.time = ztimes, misd.time = misd.times,
                       misd.rmses = misd.grmses, lsq.rmse = grmses,
                       lsq.mu.err = mu.err, misd.mu.err = norm(misd.mu.err - mu)))
# write.csv(res, "TimevN_misd.csv", row.names = F)


par(mfrow = c(1, 2))
todrop <- c(8,length(zns))
plot(zns[-todrop], misd.times[-todrop], pch = 2, 
     ylim = c(0, max(misd.times)), col='blue',
     cex.lab = 1.3,
     xlab = "Number of events", ylab = "Computation time (sec)")
lines(zns[-todrop], misd.times[-todrop], col="blue")
points(zns, ztimes, pch = 0, col="red")
lines(zns, ztimes, col="red")
legend("topleft", pch=c(0,2), lty = 1, lwd=2, col = c("red","blue"), 
       cex = 1.3, bty="n", legend = c("Least-squares", "MISD"))


plot(zns[-todrop], misd.grmses[-todrop], pch = 2, ylim = c(0, max(misd.grmses)),
     xlab = "Number of events", ylab = "RMSE", cex.lab = 1.3, col="blue")
lines(zns[-todrop], misd.grmses[-todrop], col="blue")
points(zns[-todrop], grmses[-todrop], pch = 0, col="red")
lines(zns[-todrop], grmses[-todrop], col="red")
par(mfrow = c(1, 1))
 

