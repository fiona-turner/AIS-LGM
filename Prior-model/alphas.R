alphas <- function(n,m){
  shapes <- loaddata("Prior-model/shapes")
  # find the first set of basis vectors
  A <- Xpca(shapes[[4]],3)
  # and the second set
  B <- Null(A)
  Y <- shapes[[5]]%*%B
  Y.svd <- svd(Y)
  W <- Y.svd$v[,1:3]
  A2 <- B %*% W
  V = cbind(A, A2)
  svdfour <- svd(shapes[[4]])
  # d is the square root of the eigenvalues
  # we want alpha to be able to reconstruct the forty shapes when combined with V
  # therefore we want it to be U%*%d
  # as we are using results from two separate svd, we need to construct our own U
  # U = (projected data) x (inverse of sigma)
  d <- c(svdfour$d[1:3], Y.svd$d[1:3])
  sigma <- diag(d)
  dinverse <- 1/d
  sigmainverse <- diag(dinverse)
  # project all 40 centred shapes on to the new basis
  Xproj <- shapes[[3]]%*%V
  alpha <- Xproj %*% sigmainverse %*% sigma
  return(alpha)
}
