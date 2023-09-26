#fun_2nd_L_delta() gives the second derivative of L wrt delta
#inputs:
#b: random effects of C model. usually the last element of vector z after change of variable
#delta: initial or estimated values for delta
fun_2nd_L_delta <- function(b,delta){
  lb <- length(b)
  if(lb==1){
    integrand <- 1/2 * delta^(-2) - delta^(-3) * b^2  + deriv_lnf_delta_1st(b,delta)^2
  } 
  if(lb>1){
    phi_delta <- c(1:(lb*(lb+1)/2)) 
    matrix.ind <-  matrix(NA,lb,lb)
    matrix.ind[lower.tri(matrix.ind, diag=T)] <- phi_delta
    matrix.ind <- t(matrix.ind)
    matrix.ind[lower.tri(matrix.ind, diag=T)] <-  phi_delta
    vec_delta <- stack(as.data.frame(matrix.ind))[,1]
    deriv_vecdelta_phidelta <- matrix(NA,lb^2,lb*(lb+1)/2,byrow = T)
    for(i in 1:length(phi_delta)){
      deriv_vecdelta_phidelta[,i] <- vec_delta %in% phi_delta[i]
    }
    
    second.der <- 1/2 * t(deriv_vecdelta_phidelta) %*% kronecker((invMatrix(delta)-2*invMatrix(delta)%*%b%*%t(b)%*%invMatrix(delta)),invMatrix(delta)) %*% deriv_vecdelta_phidelta
    
    integrand <- second.der + deriv_lnf_delta_1st(b,delta) %*% t(deriv_lnf_delta_1st(b,delta))

  }
  return(integrand)
}
