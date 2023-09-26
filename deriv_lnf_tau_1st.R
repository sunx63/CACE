#fun_der_1st_ln_phitau() gives the first derivative wrt tau, \partial ln(phi) / \partial phi(tau)^T for subject i in clinic j
#inputs:
#u: random effects of y model. usually the vector z of length r' after change of variable
#tau: initial or estimated values for tau

fun_der_1st_ln_phitau <- function(u,tau){
  lu <- length(u)
  if(lu==1){
    deriv_1st_lnph_tau <- 0.5 * tau^(-2) * (u^2-tau)
  }
  if(lu>1){
    #locate index of 1 and create derivative matrix of vector tau / phi tau
    phi_tau <- c(1:(lu*(lu+1)/2)) 
    matrix.ind <-  matrix(NA,lu,lu)
    matrix.ind[lower.tri(matrix.ind, diag=T)] <- phi_tau
    matrix.ind <- t(matrix.ind)
    matrix.ind[lower.tri(matrix.ind, diag=T)] <-  phi_tau
    vec_tau <- stack(as.data.frame(matrix.ind))[,1]
    deriv_vectau_phitau <- matrix(NA,lu^2,lu*(lu+1)/2,byrow = T)
    for(i in 1:length(phi_tau)){
      deriv_vectau_phitau[,i] <- vec_tau %in% phi_tau[i]
    }
    deriv_1st_lnph_tau <- 0.5 * t(deriv_vectau_phitau) %*% kronecker(invMatrix(tau),invMatrix(tau)) %*% vec(u %*% t(u)-tau)
  }
  
  return(deriv_1st_lnph_tau)
}

