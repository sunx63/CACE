
#fun_der_lnf_alpha() gives the first derivative of lnf with respect to alpha
#inputs:
#u: random effects of y model. usually the vector z of length r' after change of variable
#b: random effect of c model. usually the last element of vector z
#Alpha: alpha vector in Y model
#gamma: gamma vector in C model
#data_j: data frame of jth site
#ind.x0: column number of covariates controlled in Y model
#ind.x1: column number of covariates controlled in C model
#lambda: factor loading matrix

fun_der_lnf_alpha <- function(u,b,Alpha,gamma,data_j,ind.x0,ind.x1,lambda){
  data_j <- data_j[!is.na(data_j$Y),]
  n.j <- nrow(data_j)
  
  if(anyNA(ind.x1)){
    Xij <- as.matrix(rep(1,n.j)) #intercept only for compliance model
  } else {
    Xij <- as.matrix(cbind(rep(1,n.j),data_j[,ind.x1])) #intercept and covariates
  }
 
  Bij <- diag(rep(1,3))
  Bij <- Bij %*% lambda
  
  first.deriv <- list()
  
  for(i in 1:n.j){
    if(anyNA(ind.x0)){
      Aij <- diag(rep(1,3))
    } else {
      Aij <- cbind(diag(rep(1,3)),matrix(unlist(replicate(3,as.matrix(data_j[,ind.x0])[i,])),nrow = 3,byrow = T))
    }
    
    eta_y <- Aij %*% Alpha + Bij %*% u
    Py <- 1/(1+exp(-eta_y))
    eta_c <- as.numeric(t(as.matrix(gamma)) %*% as.matrix(Xij[i,]) + b)
     
    if(data_j$D[i]==0 & data_j$trt[i]==1){ #never taker
      first.deriv[[i]] <- (data_j$Y[i]-Py[1]) * Aij[1,]
    }

    if(data_j$D[i]==0 & data_j$trt[i]==0){ #control complier
      first.deriv[[i]] <- (exp(data_j$Y[i]*eta_y[1]-log(1+exp(eta_y[1])))*(data_j$Y[i]-Py[1])*Aij[1,] + 
                          exp(data_j$Y[i]*eta_y[2]-log(1+exp(eta_y[2]))+eta_c)*(data_j$Y[i]-Py[2])*Aij[2,]) / 
                          (exp(data_j$Y[i]*eta_y[1]-log(1+exp(eta_y[1])))+exp(data_j$Y[i]*eta_y[2]-log(1+exp(eta_y[2]))+eta_c))
    }
    if(data_j$D[i]==1 & data_j$trt[i]==1){ #trt complier
      first.deriv[[i]] <- (data_j$Y[i]-Py[3]) * Aij[3,]
    }
  }#end of loop
  
  first_der_lj_alpha <- Reduce('+', first.deriv)
 
  return(first_der_lj_alpha) 
}






