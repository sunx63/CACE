
#MtoVech() converts covariance matrix to a vector of unique elements
MtoVech <- function(CovM,d){
  phi_tau <- numeric(length = d*(d+1)/2)
  irow=0
  for(i in 1:d){
    for(j in i:d){
      irow=irow+1
      phi_tau[irow] <- CovM[i,j]
    }
  }
  return(phi_tau)
}
