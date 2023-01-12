#data load
library(ncdf4)
nc1 <- nc_open("MCMC/xkgrz.qrparm.orog.nc", write = TRUE)
controlorog <- ncvar_get(nc1, "ht")
nc_close(nc1)

#load mean of 40 shapes
mean <- read.csv('MCMC/mean.csv')
mean <- mean[,2:74]
mean <- mean[,ncol(mean):1]

o <- array(rep(0,47*96*73), dim=c(47,96,73))
design = c('1', '2', '3', '6', '8', '10', '11', '12', '13', '14', '16', '18', '19', '20', '24', '29', '31', '32', '33', '36', '37', '38', '39', '41', '42', '44', '47', '48', '49', '51', '54', '55', '56', '58', '59', '61', '64', '66', '68', '69', '70', '71', '72', '77', '78', '79', '80')

for (i in 1:length(design)){
  infile = paste("./qrparm-files/qrparm.design",design[i],".nc",sep="")
  data <- nc_open(infile)
  orog <- ncvar_get(data, "ht")
  nc_close(data)
  o[i,,] <- orog - controlorog
}

Byrdoroganom = o[,65,69]
MMoroganom = o[,61,67]
Sipleoroganom = o[,57,69]
WDCoroganom = o[,67,68]

oroganom = cbind(Byrdoroganom, MMoroganom, Sipleoroganom, WDCoroganom)

tmp = read.csv("MCMC/alpha1.csv", header = TRUE)
a1 = tmp[,2]
tmp = read.csv("MCMC/alpha2.csv", header = TRUE)
a2 = tmp[,2]
tmp = read.csv("MCMC/alpha3.csv", header = TRUE)
a3 = tmp[,2]
tmp = read.csv("MCMC/alpha4.csv", header = TRUE)
a4 = tmp[,2]
tmp = read.csv("MCMC/alpha5.csv", header = TRUE)
a5 = tmp[,2]

A = matrix(c(a1,a2,a3,a4,a5),nc=5)


nc2 = nc_open("MCMC/xkgrz.d18O_pw.climate.nc")
preind <- ncvar_get(nc2, "d18O")
nc_close(nc2)

jobs = c('xocda', 'xocdb', 'xocdc', 'xocde', 'xocdf', 'xocdg', 'xocdh', 'xocdi', 'xocdj', 'xocdk',
         'xocdl', 'xocdm', 'xocdn', 'xocdo', 'xocdp', 'xocdq', 'xocdr', 'xocds', 'xocdt', 'xocdu',
         'xocdv', 'xocdw', 'xocdx', 'xocdy', 'xocdz', 'xodja', 'xodjb', 'xodjc', 'xodjd', 'xodje',
         'xodjg', 'xodjh', 'xodji', 'xodjj', 'xodjk', 'xodjl', 'xodjm', 'xodjn', 'xodjo', 'xodjp',
         'xodjq', 'xodjr', 'xodjs', 'xodjt', 'xodju', 'xodjv', 'xodjw')

a <- array(rep(0,47*96*73), dim=c(47,96,73))

for (i in 1:length(jobs)){
  infile = paste("./HadCM3-runs/",jobs[i],".d18O_pw.climate.nc",sep="")
  data <- nc_open(infile)
  d18O <- ncvar_get(data, "d18O")
  nc_close(data)
  a[i,,] <- d18O - preind
}

#Western sites
Byrdanom = a[,65,69]
MManom = a[,61,67]
Sipleanom = a[,57,69]
WDCanom = a[,67,68]

#Eastern sites
EDCanom = a[,33,67]
EDMLanom = a[,1,67]
Fujianom = a[,11,67]
Talosanom = a[,43,66]
Tayloranom = a[,43,68]
Vostokanom = a[,29,68]

anom = cbind(Byrdanom, MManom, Sipleanom, WDCanom)
anomALL = cbind(Byrdanom, MManom, Sipleanom, WDCanom, EDCanom, EDMLanom, Fujianom, Talosanom, Tayloranom, Vostokanom)
anomEast = anomALL[,5:10]

#Western obs
Byrdobs = as.array(-7.59) #this is the value of the Byrd ice core ~21Ka 
MMobs = as.array(-4.70)
Sipleobs = as.array(-6.69)
WDCobs = as.array(-8.12)

#Eastern obs
EDCobs = as.array(-5.75)
EDMLobs = as.array(-4.68)
Fujiobs = as.array(-4.71)
Talosobs = as.array(-5.49)
Taylorobs = as.array(-1.97)
Vostokobs = as.array(-4.56)

icecoreobs = c(Byrdobs, MMobs, Sipleobs, WDCobs)
icecoreobsALL = c(Byrdobs, MMobs, Sipleobs, WDCobs, EDCobs, EDMLobs, Fujiobs, Talosobs, Taylorobs, Vostokobs)
icecoreobsEast = icecoreobsALL[5:10]
icecorenames = c("Byrd", "Mount Moulton", "Siple", "WDC", "EDC", "EDML", "Fuji", "Talos", "Taylor", "Vostok")
