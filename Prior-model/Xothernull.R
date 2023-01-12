Xothernull <- function(x,y,n){
  # x = centered dataset
  # y = existing components
  #find the null space of y, which will be the first three basis vectors here (A)
  B <- Null(y)
  #project x on to this null space. Here x is the other 36 shapes, centred
  Y <- x%*%B
  #find the svd of this projection
  Y.svd <- svd(Y)
  #the first n vectors are taken as w, the vectors to be multiplied with the nullspace
  W <- Y.svd$v[,1:n]
  #the nullspace is multiplied by these vectors to find new basis vectors for the prior model
  A2 <- B %*% W
  return(A2)
}

