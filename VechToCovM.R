#create a symmetric matrix from a vector of unique elements by row. from \phi_tau to tau

VechToCovM <- function(phi_tau,d){
  CovM <- matrix(NA,d,d)
  irow=0
  for(i in 1:d){
    for(j in i:d){
      irow=irow+1
      CovM[i,j] <- CovM[j,i] <- phi_tau[irow]
    }#end of loop j
  } #end of loop i
  return(CovM)
}#end of function


