### Compare binned OLS estimates of K against different discretizations delta(t)

source("Kestimate_functions.R")
source("simetas.R")
library(RSpectra)

######### Generate 
set.seed(375) 
bb = 250; bd = 750; 
be = 70; bf = 70


bw = 80; bm = 300; bs = 200
theta_beta = 5
a3 = 5; b3 = 2; # params for gaussian triggering
mu = 1; theta_K = 0.5; 
BigT = 1000 * 0.9 

####### Repeat estimation M times
M <- 25
Kest_rmse <- list()

ztimes <- list()
zns <- c()
Brs <- c()

xs <- seq(0, BigT, length = 10000)
delts <- seq(0.25, 5, by = .25)
Kest_avgs <- matrix(0, nrow = length(xs), ncol = length(delts))

Ws <- c()

for (rep in 1:M) {
  bb = 0.25 * BigT; bd = 0.75 * BigT;
  ba <- 1 / max(dnorm(xs, m = bb, s = be) + 0.7 * dnorm(xs, m = bd, s = bf))
  bc = 0.7 * ba
  myktest = ba * dnorm(xs, mean=bb, sd=be) + bc * dnorm(xs, mean=bd, sd=bf)
  Brs <- c(Brs, mean(myktest))
  
  if ((rep %% 5) == 0) print(rep)
  z0 = simhawk(BigT = BigT, gmi = normprod, gxy = pointxy, gt = expgt, 
               theta = list(mu = mu, K = theta_K, beta = theta_beta, b = 5,
                            ba=ba, bb=bb, bc=bc, bd=bd, be=be, bf=bf) )   
  t <- z0$t
  n <- z0$n
  zns <- c(zns, n)
  
  for (k in 1:length(delts)) {
    delt <- delts[k]
    
    ## estimate K over intervals of length I0 (continous time)
    ## which rperesents I0/delta discrete time bins
    ## In total you have T2 discrete time intervals from [0, BigT]
    ## so you have T2/I estimates of K
    I0 <- 10 
    I <- ceiling(10 / delt) 
    bin.breaks <- seq(0, BigT + delt, by = delt)
    T2 <- length(bin.breaks) - 1 
    start <- Sys.time()
    # discretize process
    W <- ceiling(T2 / I) # number of estimates of K: estimate K over 
    true.ts = bin.breaks[seq(1, T2, by = I)] # for smoothing; the times each interval corresponds to
  
    hout <- hist(z0$t, breaks = bin.breaks, plot = F)
    binCounts <- hout$counts
    
    # Construct G (from known triggering form)
    G <- matrix(0, nrow = T2 - 1, ncol = W)
    for (i in 2:T2) {
      upTo <- min(W, ceiling(i / I)) # which period the day is in
      for(j in 1:upTo) { # which periods precede bin i
        daysInJ <- seq((j-1) * I + 1, min(i - 1, j * I))
        G[i - 1, j] <- sum(dexp(i - daysInJ, 1/theta_beta) * binCounts[daysInJ])
      }
    }
    
    # Solve with ridge regularization
    sG <- svds(G, k = 1, nu = 0, nv = 0)$d # max singular value
    reg_param <- 0.01 * sG 
    L <- reg_param * diag(1, ncol = ncol(G), nrow = ncol(G))
    
    Kest0 <- (qr.solve(rbind(G, L), c(binCounts[-1] - mu * delt, rep(0, nrow(L)))))
    Kest0 <- c(Kest0, 0)
    Kest <- ksmooth(true.ts, Kest0, kernel = "normal", b = BigT/10, x.points = t)$y
    Kest <- force0(Kest)

    # Diagnostics
    myk <- ba * dnorm(t, mean=bb, sd=be) + bc * dnorm(t, mean=bd, sd=bf)
  
    Kest_avgs[,k] <- Kest_avgs[,k] +  ksmooth(t, Kest, kernel = "normal", b =  BigT/10, x.points = xs)$y
    if (rep == 1) {
      Ws <- c(Ws, W)
      ztimes[[k]] <- difftime(Sys.time(), start, units = "secs")
      Kest_rmse[[k]] <- mse(myk - Kest)
    } else {
      ztimes[[k]] <- c(ztimes[[k]], difftime(Sys.time(), start, units = "secs"))
      Kest_rmse[[k]] <- c(Kest_rmse[[k]], mse(myk - Kest))
    }
  }
}


## plot
par(mfrow = c(1, 2))
plot(delts, unlist(lapply(ztimes, mean)), 
     xlab = expression(Delta), ylab = "Time (secs)", pch = 3)
plot(delts, unlist(lapply(Kest_rmse, mean)), 
     xlab = expression(Delta), ylab = "RMSE", pch = 2)
par(mfrow = c(1, 1))



# Plot estimates
ndt <- length(delts)
colfunc <- colorRampPalette(c("orange", "red"))
cols <- colfunc(length(delts))
myk <- ba * dnorm(xs, mean=bb, sd=be) + bc * dnorm(xs, mean=bd, sd=bf)
plot(xs, myk, xlab = "w", ylab = "K(w)", t = "l",
     ylim = range(myk),
     main = paste0("K(t) ~ N(",bb,",",be,") + N(",bd,",",bf, ")"))
for (i in 1:length(delts)) {
  Ksm <- ksmooth(xs, Kest_avgs[,i] / M, kernel = "normal", b = bw.nrd(xs), x.points = xs)$y
  lines(xs, Ksm, col = cols[i], lty = 2)
}

legend("topright", col = cols[c(1, 7, ndt)], lty = 2, 
       legend = paste("I =", c(1,7,ndt)))

