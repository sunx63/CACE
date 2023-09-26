
#AGHQ() approximate integrals with a sum of weighted Q abscissas
#inputs: 
#Q: number of abscissas
#Wphi: weights * probabilities of a vector of z (abscissas after change of variables) that follow multivariate standard normals distribution
#f: g()f(Y,C), g() is the function to be approximated, f(Y,C) is the joint distribution of outcome Y and compliance 
#z: a vector or a matrix of abscissas after change of variable
#r: reduced dimension of random effects of Y model
#k: reduced dimension of random effects of C model. it is 1 for one-sided noncompliance
#Alpha: alpha vector in Y model
#gamma: gamma vector in C model
#data_j: data frame of jth site
#ind.x0: column number of covariates controlled in Y model
#ind.x1: column number of covariates controlled in C model
#lambda: factor loading matrix
#tau: variance (matrix) of random effects of Y model
#delta: varianceof random effect of C model
#M: derivative of factor loading matrix

#Output: approximation of an integral

AGHQ <- function(Q,Wphi,f,z,r,k,Alpha,gamma,data_j,ind.x0,ind.x1,lambda,tau,delta,M){ #f:function. add z as an additional arguement
  p <- Q^(r+k)
  L=0
  if(missing(M)){
    if(missing(tau) & missing(delta)){
      for(q in 1:p){ 
        L <- L + f(z[q,1:r],z[q,(r+1):(r+k)],Alpha,gamma,data_j,ind.x0,ind.x1,lambda)*Wphi[q] 
      }
    }
    if (missing(tau) & !missing(delta)) {
      for(q in 1:p){ 
        L <- L + f(z[q,1:r],z[q,(r+1):(r+k)],Alpha,gamma,data_j,ind.x0,ind.x1,lambda,delta)*Wphi[q] 
      }
    }
    if (!missing(tau) & missing(delta)) {
      for(q in 1:p){ 
        L <- L + f(z[q,1:r],z[q,(r+1):(r+k)],Alpha,gamma,data_j,ind.x0,ind.x1,lambda,tau)*Wphi[q] 
      }
    }
    if (!missing(tau) & !missing(delta)) {
      for(q in 1:p){ 
        L <- L + f(z[q,1:r],z[q,(r+1):(r+k)],Alpha,gamma,data_j,ind.x0,ind.x1,lambda,tau,delta)*Wphi[q] 
      }
    }
  } else {#end of missing(M)
    if (missing(tau) & missing(delta)) {
      for(q in 1:p){ 
        L <- L + f(z[q,1:r],z[q,(r+1):(r+k)],Alpha,gamma,data_j,ind.x0,ind.x1,lambda,M)*Wphi[q] 
      }
    }
    if (!missing(tau) & missing(delta)) {
      for(q in 1:p){ 
        L <- L + f(z[q,1:r],z[q,(r+1):(r+k)],Alpha,gamma,data_j,ind.x0,ind.x1,lambda,tau,M)*Wphi[q] 
      }
    }
    if (missing(tau) & !missing(delta)) {
      for(q in 1:p){ 
        L <- L + f(z[q,1:r],z[q,(r+1):(r+k)],Alpha,gamma,data_j,ind.x0,ind.x1,lambda,delta,M)*Wphi[q] 
      }
    }
  } #end of !missing(M)
  return(L)
}


