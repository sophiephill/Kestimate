## General ETAS simulator where the user can input densities.

simhawk = function(x1=1, y1=1, T=100, rho=unifrho, gt=powergt, gxy=powerxy, 
                   gmi=expprod, mdensity=expmag, sor=1, keep=1, verbose = FALSE,
                   theta = list(mu=8.0,K=.3,c=.2,p=2.3,a=.1,d=.3,q=2.7,b=.7)){
  ##### THIS IS FOR SIMULATING A HAWKES PROCESS WITH 
  ##### lambda(t,x,y) = mu rho(x,y) + 
  ##### SUM gmi(m_i) gt(t-t_i) gxy(x-xi,y-yi; mi),
  ##### on a space S = [0,x1] x [0,y1] (km), in time [0,T], 
  ##### background temporal rate mu and spatial density rho(x,y),
  ##### triggering density gt(t-t_i) gxy(x-xi, y-yi; mi),
  ##### productivity gmi(m_i),
  ##### and magnitude density mdensity(m).
  ##### sor = 1 outputs the points in chronological order.
  ##### keep = 1 means only keep the ones within the space time window.
  ##### Both gt and gxy must be densities, so that if mu = 1/(x1y1),
  ##### then the integral of lambda over the space time region = mu T + SUM gmi(m_i,m0).
  ##### Thus the ETAS parameter K is included in gmi.
  ##### If no magnitudes are desired, just let gmi = K. 
  y = bgpts(x1,y1,T, rho, mdensity,theta) ## lay down the background points.
  bgt = y$t
  if(verbose) cat(y$n,"mainshocks.\n") 
  calcbr = 0
  #calcbr = mean(gmi(mdensity(1000000,theta),theta)) ## calculate branching ratio. Stop if br > 1. 
  if (verbose) cat("branching ratio is ", calcbr,"\n")
  if(calcbr>1.0){
    if (verbose)  cat("error, branching ratio = ", calcbr, " > 1.")
    return(0)
  }
  stop1 = 0
  if(y$n < 0.5) stop1 = 2
  if (verbose) cat("aftershocks by generation\n")
  w = y
  while(stop1 < 1) {
    z = aft(w,x1,y1,T,gt,gxy,gmi,mdensity,theta) ## place aftershocks down around y.
    if (verbose) cat(z$n," ")
    if(z$n > 0.5){
      y = combine1(y,z)
      w = z
      if(min(z$t) > T) stop1 = 2
    }
    if(z$n < 0.5) stop1 = 1
  }
  y$Kn <- c(y$Kn, z$Kn)
  if(keep==1) y = keep1(y,x1,y1,T) ## to keep just the pts in the window.
  if(sor==1) y = sort1(y) ## to have the points sorted chronologically.
  y
}

## br = INT gmi(m) mdensity(m) dm, from m = m0 to infinity.

bgpts = function(x1,y1,T, rho, mdensity,theta){
  z1 = list()
  n = rpois(1,theta$mu*T)
  z1$n = n
  xy = rho(n,x1,y1)
  z1$lon = xy[,1]
  z1$lat = xy[,2]
  z1$t = sort(runif(n)*T)
  z1$m = mdensity(n,theta)
  z1$ztimes = c()
  z1$Kn = c()
  z1
}


# need to modify aft() so it accepts time dependent gmi (productivity)
aft = function(y,x1,y1,T,gt,gxy,gmi,mdensity,theta){
  ## place aftershocks around y.
  z1 = list()
  n2 = gmi(y$m,y$t,theta) ## vector of number of aftershocks for each mainshock.
  for(i in 1:length(n2)){
    if(n2[i] > 0.5){
      b1 = gt(n2[i],theta)
      z1$ztimes = c(z1$ztimes, b1)
      z1$t = c(z1$t, b1 + y$t[i])
      xy = gxy(n2[i], y$m[i])
      z1$lon = c(z1$lon, xy[,1] + y$lon[i])
      z1$lat = c(z1$lat, xy[,2] + y$lat[i])
      z1$m = c(z1$m, mdensity(n2[i],theta))
    }
  }
  z1$Kn = n2
  z1$n = sum(n2)
  return(z1)
}

combine1 = function(y,z){
  z1 = list()
  z1$t = c(y$t,z$t)
  z1$n = y$n + z$n
  z1$m = c(y$m,z$m)
  z1$lat = c(y$lat,z$lat)
  z1$lon = c(y$lon,z$lon)
  z1$Kn = c(y$Kn, z$Kn)
  z1$ztimes = c(y$ztimes, z$ztimes)
  z1
}

keep1 = function(y,x1,y1,T){
  ## keep only the pts of y that are within the space time window [0,x1] x [0,y1] x [0,T].
  keeps = c(1:length(y$t))[(y$t<T)&(y$lon<x1)&(y$lat<y1)&(y$lon>0)&(y$lat>0)]
  y$t = y$t[keeps]
  y$m = y$m[keeps]
  y$lon = y$lon[keeps]
  y$lat = y$lat[keeps]
  y$n = length(keeps)
  y$Kn = y$Kn[keeps]
  y
}

sort1 = function(y){
  ## sort the pts chronologically.
  ord2 = order(y$t)
  y$t = y$t[ord2]
  y$m = y$m[ord2]
  y$lon = y$lon[ord2]
  y$lat = y$lat[ord2]
  y$Kn = y$Kn[ord2]
  y
}

## rho takes an integer n and x1 and y1 and outputs a matrix of n locations of mainshocks.
## gt takes an integer n and outputs a vector of n nonnegative times since mainshock.
## gxy takes an integer n and magnitude m and outputs a matrix of n locs from mainshock.
## gmi takes a vector of mags m and outputs a vector of number of aftershocks per mainshock.   
## mdensity takes an integer n and lower mag threshold m0 and outputs a vector of n magnitudes.

## Below are examples of functions rho, gt, gxy, gmi, and mdensity.

unifrho = function(n,x1,y1){
  ## Uniform spatial background density rho on [0,x1] x [0,y1].
  x = runif(n,min=0,max=x1)
  y = runif(n,min=0,max=x1)
  cbind(x,y)
}

expmag = function(n,theta, m0=3.5){
  ## exponential magnitude density mdensity with minimum m0 and mean m0+b1.
  rexp(n,rate=1/theta$b) + m0
}

powergt = function(n,theta){
  ## power law triggering function in time gt
  ## f(u) = (p-1) c^(p-1) (u+c)^-p.
  v = runif(n)
  theta$c*(1-v)^(1/(1-theta$p)) - theta$c
}

## Notes for powergt.
## if v = runif(1), then new time is found by letting v = F(t) and solving for t.
## F(t) = INT from 0 to t of f(u) du = (p-1) c^(p-1) (u+c)^(1-p) / (1-p)
## from u = 0 to t
## = -c^(p-1) (t+c)^(1-p) + c^(p-1) c^(1-p) = 1 - c^(p-1) (t+c)^(1-p).
## Setting v = 1 - c^(p-1) (t+c)^(1-p) and solving for t, we get
## c^(p-1) (t+c)^(1-p) = 1-v.
## (t+c)^(1-p) = (1-v) c^(1-p).
## t+c = c (1-v)^{1/(1-p)}.
## t = c (1-v)^{1/(1-p)} - c.

expgt = function(n,theta){
  ## exponential triggering function in time gt, with mean beta2. 
  ## f(u) = beta e^(-betau).
  rexp(n,rate=1/theta$beta)
}

normg = function(n, theta){
  ## normal triggering function in time gt
  abs(rnorm(n,mean=theta$a3,sd=theta$b3))
}


powerxy = function(n,m,theta){
  ## power law triggering in space according to ETAS (2.3), gxy, of
  ## Ogata (1998). See http://wildfire.stat.ucla.edu/pdflibrary/ogata98.pdf . 
  ## Here the density does not depend on magnitude of the mainshock. 
  ## ∫f(x,y)dxdy = 1 = ∫h(r)rdrdø = 2π∫h(r)rdr. 
  ## h(r) = c (r^2 + d)^(-q). 
  ## ∫ h(r)rdr = c(r^2+d)^(1-q)/(2-2q),r=0to∞. For q > 1, this is 0+cd^(1-q)/(2q-2).
  ## So c = (q-1)d^(q-1)/π. 
  v = runif(n)
  dist1 = sqrt(theta$d*(1-v)^(1/(1-theta$q))-theta$d)
  thet1 = runif(n)*2*pi
  x = cos(thet1)*dist1
  y = sin(thet1)*dist1
  cbind(x,y)
}

expxy = function(n,m,theta){
  ## exponential triggering in space. f(r) = alpha/pi exp(-alpha r^2). 
  ## Here the density does not depend on magnitude of the mainshock. 
  ## To see that this is a density, 
  ## ∫f(x,y)dxdy = ∫f(r)rdrdø = 2π∫f(r)rdr 
  ## = 2alpha ∫ exp(-alpha r^2) r dr = -exp(-alpha r^2) ,r=0to∞, = 0+1, for alpha>0.  
  v = rexp(n,rate=theta$alpha)
  dist1 = sqrt(v)
  thet1 = runif(n)*2*pi
  x = cos(thet1)*dist1
  y = sin(thet1)*dist1
  cbind(x,y)
}

## gaussian renewal productivity 
normprod = function(m, t, theta) {
  pmean <- theta$ba*dnorm(t,mean=theta$bb,sd=theta$be) + theta$bc*dnorm(t,mean=theta$bd,sd=theta$bf)
  rpois(length(t), pmean)
}

normprod2 = function(m, t, theta) {
  pmean <- theta$bw * dnorm(t, m = theta$bm, s = theta$bs)
  rpois(length(t), pmean)
}

pointprod = function(m, t, theta) rpois(length(t), theta$K) ## Here each point has productivity K. 

## Examples. 
# b = simhawk()
# b = simhawk(T=10)
# b = simhawk(T=10,theta = list(mu=8.0,K=.3,c=.2,p=2.3,a=.1,d=.3,q=2.8,b=.7))
# mags = (b$m - min(b$m))/diff(range(b$m))
# par(mfrow=c(1,2))
# plot(b$lon,b$lat,pch=3,cex = 1+3*mags,xlab="longitude (km)",ylab="latitude (km)")
# plot(b$t,b$m,pch=3,xlab="time (days)", ylab = "magnitude")
# 
# # install.packages("gplm")
# library(gplm)
# library(MASS)
# t2 = powergt(100000,theta=list(c=.2,p=2.3))
# hist(t2[t2<1],nclass=100,prob=T,main="",xlab="time between pts")
# lines(kde(t2[t2<1]),col="green")
# 
# g2 = powerxy(100000,1,theta=list(d=.3,q=2.8))
# g3 = kde2d(g2[,1],g2[,2],lims=c(-1,1,-1,1))
# image(g3)
# contour(g3,add=T)
