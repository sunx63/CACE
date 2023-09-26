#\partial^2 L / \partial tau^2

#fun_2nd_L_tau() gives the second derivative of L wrt tau
#inputs:
#u: random effects of y model. usually the vector z of length r' after change of variable
#tau: initial or estimated values for tau
fun_2nd_L_tau <- function(u,tau){
    lu <- length(u)
    if(lu==1){
     integrand <- 1/2 * tau^(-2) - tau^(-3) * u^2  + fun_der_1st_ln_phitau(u,tau)^2
    } 
    if(lu>1){
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
      
      deriv_vectau_phitau <- ifelse(deriv_vectau_phitau==T, 1, 0) #newly added
      
      product_u_tau <- invMatrix(tau)%*%u%*%t(u)%*%invMatrix(tau) #newly added for solving asymmertric Htau

      second.der <- 1/2 * t(deriv_vectau_phitau) %*% kronecker((invMatrix(tau)-2*product_u_tau),invMatrix(tau)) %*% deriv_vectau_phitau
      
      integrand <- second.der + fun_der_1st_ln_phitau(u,tau) %*% t(fun_der_1st_ln_phitau(u,tau))
     # integrand <- (integrand + t(integrand))/2 #fix asymmetric Ht issue
    }
   return(integrand)
}

