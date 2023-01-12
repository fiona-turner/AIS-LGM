Xpca <- function(x,n){
  #perform svd on the matrix
  svdfour <- svd(x)
  #take the first n singular column vectors
  A <- svdfour$v[,1:n]
  #these are our first n basis vectors
  return(A)
}
