source("MCMC/Data-load.R")
library(rstan)


toy_dataALL <- 
  list(N=dim(anomALL)[[1]], #47 - number of runs
       M=dim(anomALL)[[2]],#10 - number of sites
       y1=anomALL, # 47 * 10 - HadCM3 output at 10 sites for 47 runs
       x1 = A, #47*5
       yobs=icecoreobsALL, #10*1
       xx1 = A
  )

#run emulator to recover a HadCM3 simulation
toy_dataALLtest <- 
  list(N=dim(anomALL[2:47,])[[1]], 
       M=dim(anomALL[2:47,])[[2]],
       y1=anomALL[2:47,], 
       x1 = A[2:47,],
       yobs=anomALL[1,],
       xx1 = A[2:47,]
  )

fitALL <- stan(file="MCMC/stancodewithgpmultiALL.stan", data = toy_dataALL,  iter = 5000, chains = 4)
chainsALL = extract(fitALL)

fitALLtest <- stan(file="MCMC/stancodewithgpmultiALL.stan", data = toy_dataALLtest,  iter = 5000, chains = 4)
chainsALLtest = extract(fitALLtest)

hist(chainsALL$alpha)

par(mfrow=c(3,4))
for(i in 1:10){
  plot(density(chainsALL$ypred[,i]), main=icecorenames[i], xlab=expression(paste(delta^{18}, "O estimation")), xlim=c(-12,2), ylim=c(0,0.5), ylab='', cex.main=2, cex.lab=2, cex.axis=2, lwd=2)
  abline(v=icecoreobsALL[i], lwd=3, col=2)
}

par(mfrow=c(3,4))
for(i in 1:10){
  plot(density(chainsALLtest$ypred[,i]), xlab=expression(paste(delta^{18}, "O estimation")), main=icecorenames[i], xlim=c(-16,8), ylim=c(0,0.5), ylab='', cex.main=2, cex.lab=2, cex.axis=2, lwd=2)
  abline(v = anomALL[47,i], col='red', lwd = 3)
}

traceplot(fitALL, pars='lp__')
traceplot(fitALL, pars='xobs')
traceplot(fitALL, pars='ypred')

par(mfrow=c(2,3))
hist(chainsALL$xobs[,1])
hist(chainsALL$xobs[,2])
hist(chainsALL$xobs[,3])
hist(chainsALL$xobs[,4])
hist(chainsALL$xobs[,5])

priorsd <- c(0.5,0.5,0.6,0.5,1)
par(mfrow=c(2,3))
for(i in 1:5) {
  hist(rnorm(1000, 0, priorsd[i]), probability=TRUE, xlab=expression(paste(theta, i)), main="", xlim=c(-2,2), ylim=c(0,2.5), ylab='', cex.lab=2, cex.axis=2)
  lines(density(chainsALLtest$xobs[,i]),col='red')
  abline(v=A[1,i], lwd=2, col=2)
}

par(mfrow=c(2,3))
hist(chainsALLtest$xobs[,1], probability=TRUE, xlab=expression(theta[1]), main="", xlim=c(-2,2), ylim=c(0,2.5), ylab='', cex.lab=2, cex.axis=2)
f<- function(x) dnorm(x, 0, priorsd[1])
curve(f, add=TRUE, col='red')
abline(v=A[1,1], lwd=2, col=2)
hist(chainsALLtest$xobs[,2], probability=TRUE, xlab=expression(theta[2]), main="", xlim=c(-2,2), ylim=c(0,2.5), ylab='', cex.lab=2, cex.axis=2)
f<- function(x) dnorm(x, 0, priorsd[2])
curve(f, add=TRUE, col='red')
abline(v=A[1,2], lwd=2, col=2)
hist(chainsALLtest$xobs[,3], probability=TRUE, xlab=expression(theta[3]), main="", xlim=c(-2,2), ylim=c(0,2.5), ylab='', cex.lab=2, cex.axis=2)
f<- function(x) dnorm(x, 0, priorsd[3])
curve(f, add=TRUE, col='red')
abline(v=A[1,3], lwd=2, col=2)
hist(chainsALLtest$xobs[,4], probability=TRUE, xlab=expression(theta[4]), main="", xlim=c(-2,2), ylim=c(0,2.5), ylab='', cex.lab=2, cex.axis=2)
f<- function(x) dnorm(x, 0, priorsd[4])
curve(f, add=TRUE, col='red')
abline(v=A[1,4], lwd=2, col=2)
hist(chainsALLtest$xobs[,5], probability=TRUE, xlab=expression(theta[5]), main="", xlim=c(-2,2), ylim=c(0,2.5), ylab='', cex.lab=2, cex.axis=2)
f<- function(x) dnorm(x, 0, priorsd[5])
curve(f, add=TRUE, col='red')
abline(v=A[1,5], lwd=2, col=2)

par(mfrow=c(3,4))
for(i in 1:10){
  hist(chainsALL$ypred[,i], main=icecorenames[i], xlab="MCMC isotope prediction from all sites")
  abline(v=icecoreobsALL[i], lwd=2, col=2)
}


pairs(chainsALL$xobs)
cor(chainsALL$xobs)
pairs(cbind(chainsALL$xobs, chainsALL$beta1))

library("invgamma")
par(mfrow=c(2,3))
for(i in 1:5){
  hist(chainsALL$rho[,i], probability=TRUE, ylim=c(0,1.55), main="", xlab=expression(rho))
  f<- function(x) dinvgamma(x, 5, 5)
  curve(f, add=TRUE, col='red')
}

det(cov(chainsALL$xobs))
#[1] 1.594461e-07
det(cov(A))
#[1] 1.863949e-05

par(mfrow=c(3,4))
for(i in 1:10){
  plot(density(chainsALL$sigma1[,i]), xlab="Emulator variance", main=icecorenames[i], xlim=c(0.4,2.2))
}
