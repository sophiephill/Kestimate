### As a comparison to MISD, used in Kcompare paper

source("simetasmay2017.r")
source("misd.R")
source("Kestimate_functions.R")
library(RSpectra)

#### generate 
X1 = 1; Y1 = 1
a3 = 7; b3 = 2 # gaussian triggering parameters
theta_beta = 5; # exp triggering parameters
mu = 1
theta_K = 0.5 ## productivity
T = 1000

#### Binning parameters
t.delta <- 0.5
day.breaks <- seq(0, T + t.delta, by = t.delta)
ndays <- length(day.breaks) - 1
W2.true <- 20 # the cont. interval you want to measure g on
W2 <- W2.true / t.delta # on many discrete time intervals this represents
true.ts <- day.breaks[1:W2 + 1] # the continuous times each interval represents


#### Repeat estimation procedure M times
M <- 25

grmses <- double(M)
mu.err <- double(M)
Kemp.err <- double(M)
mu.ests <- double(M)

rparams <- double(M)
zns <- double(M)
ztimes <- double(M)

gavg <- double(W2)

misd.nbins <- 15
grange <- seq(0, W2.true, length = misd.nbins + 1)
delt <- diff(grange)[1]

misd.gavg <- double(misd.nbins)
misd.grmses <- double(M)
misd.mu.ests <- double(M)
misd.times <- double(M)

set.seed(404) 
for (rep in 1:M) {
  print(rep)
  z0 = simhawk(T = T, gmi = pointprod, gxy = pointxy, gt = expgt, 
               theta = list(mu = mu, 
                            b = 5, # for mdensity
                            K = theta_K, # for uniform productivity
                            beta = theta_beta, # exp triggering
                            a3 = a3, b3 = b3, # gaussian triggering
                            ba=ba, bb=bb, bc=bc, bd=bd, be=be, bf=bf
               ) )
  n <- z0$n
  zns[rep] <- n
  bout <- hist(z0$t, day.breaks, plot = FALSE) # bin unique event times into "days"
  dailyCounts <- bout$counts
  
  start <- Sys.time()
  A <- matrix(0, nrow = ndays - 1, ncol = W2)
  for (i in 2:ndays) { 
    upTo <- min(i - 1, W2) # how many days before
    for(j in 1:upTo) { 
      dayInd <- i - j
      weekInd <- ceiling(dayInd / I) # which period the day is in
      A[i - 1, j] <-  theta_K * dailyCounts[dayInd]
    }
  }
  A <- cbind(1, (A))
  
  sA <- 0.1 * svds(A, k = 1, nu = 0, nv = 0)$d # max sing val of A
  rparams[rep] <- sA
  
  L <- diag(c(0, rep(1, W2)), ncol = W2 + 1) # ridge
  gest0 <- qr.solve(rbind(A, sA * L), c(dailyCounts[-1], rep(0, nrow(L))) )
  
  mu.ests[rep] <- gest0[1]
  
  gest <- ksmooth(true.ts, gest0[-1], b = bw.nrd(true.ts), kernel = "normal", x.points = true.ts)$y
  if (all(gest <= 0)) stop("all negative")
  gest <- (gest) / sum(t.delta * pmax(gest, 0))
  gavg <- gavg + gest
  
  gtrue <- dnorm(true.ts, m = a3, s = b3) / sum(t.delta * dnorm(true.ts, m = a3, s = b3))
  # gtrue <- dexp(true.ts, 1/theta_beta) / sum(t.delta * dexp(true.ts, 1/theta_beta))
  grmses[rep] <- norm(gtrue - gest)
  
  ztimes[rep] <- difftime(Sys.time(), start, units = "secs")
  
  #### Compare to misd
  misd.start <- Sys.time()
  misd.out <- misd(z0$t, grange = grange, num_iter = 1e3, tot_time = T, verbose = F)
  misd.gest <- misd.out$g / sum(delt * misd.out$g)

  misd.gtrue <- dnorm(grange, m = a3, s = b3) / sum(delt * dnorm(grange, m = a3, s = b3))
  # misd.gtrue <- dexp(grange[-1], 1/theta_beta) / sum(delt * dexp(grange[-1], 1/theta_beta))

  misd.gavg <- misd.gavg + misd.gest
  misd.grmses[rep] <- norm(misd.gtrue - misd.gest)
  misd.mu.ests[rep] <- misd.out$mu
  misd.times[rep] <- difftime(Sys.time(), misd.start, units = "secs")
}

sqrt(mean((mu - mu.ests/t.delta)^2))


#### Table data
data.frame(list(density = paste0("Norm(", a3, ",", b3, ")"), 
                #density = paste0("g(t) ~ Exp(", theta_beta, ")"),
                gavgRmse = norm(gtrue - gavg / sum(t.delta * pmax(gavg, 0))),
                avgMisdRmse = norm(misd.gtrue - misd.gavg / sum(delt * misd.gavg)),
                MeangRMSE = mean(grmses), MeanMisdRmse = mean(misd.grmses),
                MeanMuRMSE = sqrt(mean((mu.ests - mu)^2)),
                MeanMisdMuRMSE = sqrt(mean((misd.mu.ests - mu)^2)),
                MeanTime = mean(ztimes), MeanMisdTime = mean(misd.times),
                MeanN = mean(zns)))


#### Plot
gtrue <- dexp(true.ts, 1/theta_beta) / sum(t.delta * dexp(true.ts, 1/theta_beta))

plot(true.ts, gtrue, type = "l", col = 1, lwd = 2,
     ylim = c(0, max(gtrue, (gavg) / sum(t.delta * pmax(gavg, 0)),
                      (misd.gavg) / sum(delt * misd.gavg))),
     xlab = "lag", ylab = "g(t)", cex.lab = 1.25)
lines(true.ts, (gavg) / sum(t.delta * pmax(gavg, 0)), lty = 2, lwd = 2, col = 2)
lines(grange[-misd.nbins], (misd.gavg) / sum(delt * misd.gavg), lty = 2, lwd = 2, col = "blue")

legend("topright", lty = c(1, 2, 2), lwd=2, col = c(1, "red","blue"), 
       cex = 1.3,
       bty="n", legend = c("True", "Least-squares", "MISD"))

