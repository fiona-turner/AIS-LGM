# analysis using functions
library("MASS")
source("Prior-model/loaddata.R")
source("Prior-model/Xpca.R")
source("Prior-model/Xothernull.R")
source("Prior-model/reconstruction.R")
source("Prior-model/alphas.R")

#load the data
shapes <- loaddata("Prior-model/shapes")
# this gives output in a list: X,mu,Xc,Xfourc,Xotherc


### method 1
# perform PCA on all 40 shapes
V1 <- Xpca(shapes[[3]],6)
# reconstruct the data using this new basis
Xrecon1 <- reconstruct(shapes[[3]],V1)
# find the errors
rmse_40 <- reconstruction_error(shapes[[3]],Xrecon1)
rmse_40

### method 2
# perform PCA on the first 4, and then the other 36, shapes
PC1 <- Xpca(shapes[[4]],3)
PC2 <- Xpca(shapes[[5]],3)
PC <- cbind(PC1, PC2)
# orthogonalise these components
V2 <- qr.Q(qr(PC))
# project the data on to the new basis
Xrecon2 <- reconstruct(shapes[[3]],V2)
# find the errors
rmse_33 <- reconstruction_error(shapes[[3]],Xrecon2)
rmse_33

### method 3
#find the first n basis vectors using X4 centred
A <- Xpca(shapes[[4]],3)
#now find the further m basis vectors, using X36 centred
A2 <- Xothernull(shapes[[5]],A,2)
#bind them in to a basis matrix
V3 = cbind(A, A2)
# project the data on to the new basis
Xrecon3 <- reconstruct(shapes[[3]],V3)
# find the errors
rmse_null <- reconstruction_error(shapes[[3]],Xrecon3)
rmse_null


# compare the three methods
rmse <- cbind(rmse_40,rmse_33,rmse_null)
rmse
rmse_means <- colMeans(rmse)
#rmse_40   rmse_33 rmse_null 
#13.52467  14.14710  14.04517  
# compare for Xfour
rmse_topfour <- colMeans(rmse[1:4,])
#rmse_40   rmse_33 rmse_null
#31.28293  21.54119  21.83490


#reconstruct first 4 shapes with A
Xreconfour <- project(shapes[[4]],A)
rmse_four <- reconstruction_error(shapes[[4]],Xreconfour)
#why are these values not close to 0?
#and not the same as prcomp
#height5g.csv  height6g.csv heightpip.csv heighttar.csv 
#1.585066      6.760732     39.661273     49.210842 
#also,why are the rmse for prcomp not zero? what's happening?
pca <- prcomp(Xfourc)
Xfourcc <- sweep(Xfourc,2,colMeans(Xfourc))
svdfourc <- svd(Xfourcc)
svdfourc$d/pca$sdev
#prcomp is recentring the data before performing pca
pca2 <- prcomp(shapes[[4]], center=FALSE)
A/pca2$rotation[,1:3]
# test the accuracy of A in reconstructing X4
Xrecon4 <- project(shapes[[4]],A)
rmse_4 <- reconstruction_error(shapes[[4]],Xrecon4)
#height5g.csv  height6g.csv heightpip.csv heighttar.csv 
#1.585066      6.760732     39.661273     49.210842

############################################################################
## trivial example
Z <- matrix(rnorm(9),nrow=3)
Zbar = colMeans(Z)
# centre it
Zc = sweep(Z,2, Zbar)
#perform svd on the matrix
svdz <- svd(Zc[1:2,])
#take the first n singular column vectors
A <- svdz$v
B <- Null(A)
Y <- Zc[3,]%*%B
#find the svd of this projection
Y.svd <- svd(Y)
#the first n vectors are taken as w, the vectors to be multiplied with the nullspace
W <- Y.svd$v
#the nullspace is multiplied by these vectors to find new basis vectors for the prior model
A2 <- B %*% W
V <- cbind(A,A2)
Zrecon = Zc%*%V%*%t(V)
rmse_Z <- sqrt(rowMeans((Zc-Zrecon)^2))

d <- c(svdz$d, Y.svd$d[1:2])
sigma <- diag(d)
dinverse <- 1/d
sigmainverse <- diag(dinverse)
# project all 40 centred shapes on to the new basis
Zproj <- Z%*%V
alpha <- Zproj %*% sigmainverse %*% sigma

