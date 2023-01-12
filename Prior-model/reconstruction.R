#project the data on to the new basis
#X is the original data, basis is the new basis vectors 
project<- function(X,Y){
  Xproj <- X%*%Y
  return(Xproj)
}

reconstruct <- function(X, basis){
  Xrecon = X%*%basis%*%t(basis)
  return(Xrecon)
}

#find the errors in the projection
#X is the original data, Y is the reconstructed data
reconstruction_error <- function(X, Y){
  rmse <- sqrt(rowMeans((X-Y)^2))
  return(rmse)
}


# 
# Prepare data by centering
# Create basis
# e.g. a) Plain pca
#      b) PCA on two separate datasets
#      c) Fancy method

# Test reconstruction errors
