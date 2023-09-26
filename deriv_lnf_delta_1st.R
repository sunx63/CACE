#deriv_lnf_delta_1st() gives the first derivative of lnf wrt delta
#inputs:
#b: random effects of C model. usually the last element of vector z after change of variable
#delta: initial or estimated values for delta

deriv_lnf_delta_1st <- function(b,delta){
  lb <- length(b)
  if(lb==1){
    deriv_1st_lnph_delta <- 0.5 * delta^(-2) * (b^2-delta)
  }
  if(lb>1){
    #locate index of 1 and create derivative matrix of vector delta / phi delta
    phi_delta <- c(1:(k*(k+1)/2)) 
    matrix.ind <-  matrix(NA,k,k)
    matrix.ind[lower.tri(matrix.ind, diag=T)] <- phi_delta
    matrix.ind <- t(matrix.ind)
    matrix.ind[lower.tri(matrix.ind, diag=T)] <-  phi_delta
    vec_delta <- stack(as.data.frame(matrix.ind))[,1]
    deriv_vecdelta_phidelta <- matrix(NA,k^2,k*(k+1)/2,byrow = T)
    for(i in 1:length(phi_delta)){
      deriv_vecdelta_phidelta[,i] <- vec_delta %in% phi_delta[i]
    }
    deriv_1st_lnph_delta <- 0.5 * t(deriv_vecdelta_phidelta) %*% kronecker(invMatrix(delta),invMatrix(delta)) %*% vec(b %*% t(b)-delta)
  }
  
  return(deriv_1st_lnph_delta)
}
