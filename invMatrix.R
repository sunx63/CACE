#2022-7-21
#find the inverse of a matrix
invMatrix <- function(matrix){
  L <- chol(matrix)
  invM <- solve(L)%*%t(solve(L))
  return(invM)
}
