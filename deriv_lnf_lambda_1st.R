
#fun_der_lnf_lambda() gives the first derivative of lambda \pratial lnf(Yj,c|u,b) / \partial lambda for subject i in clinic j
#inputs:
#u: random effects of y model. usually the vector z of length r' after change of variable
#b: random effect of c model. usually the last element of vector z
#Alpha: alpha vector in Y model
#gamma: gamma vector in C model
#data_j: data frame of jth site
#ind.x0: column number of covariates controlled in Y model
#ind.x1: column number of covariates controlled in C model
#lambda: factor loading matrix
#M is the direvative of factor loading matrix

fun_der_lnf_lambda <- function(u,b,Alpha,gamma,data_j,ind.x0,ind.x1,lambda,M){
  data_j <- data_j[!is.na(data_j$Y),]
  n.j <- nrow(data_j)
  
  if(anyNA(ind.x1)){
    Xij <- as.matrix(rep(1,n.j)) #intercept only for compliance model
  } else {
    Xij <- as.matrix(cbind(rep(1,n.j),data_j[,ind.x1])) #intercept and covariates
  }

  Bij <- diag(rep(1,3))

  first.deriv <- list()
  
  for(i in 1:nrow(data_j)){
    if(anyNA(ind.x0)){
      Aij <- diag(rep(1,3))
    } else {
      Aij <- cbind(diag(rep(1,3)),matrix(unlist(replicate(3,as.matrix(data_j[,ind.x0])[i,])),nrow = 3,byrow = T))
    }
    
    eta_y <- Aij %*% Alpha + Bij %*% lambda %*% u
    Py <- 1/(1+exp(-eta_y))
    eta_c <- as.numeric(t(as.matrix(gamma)) %*% as.matrix(Xij[i,]) + b)
    
    deriv.etay.lambda <- t(Bij) %*% M %*% u
    

    if(data_j$D[i]==0 & data_j$trt[i]==1){ #never taker
     first.deriv[[i]] <- (data_j$Y[i]-Py[1]) * deriv.etay.lambda[1]
    }
    
    if(data_j$D[i]==0 & data_j$trt[i]==0){ #control complier
     D00 <- exp(data_j$Y[i]*eta_y[1]-log(1+exp(eta_y[1]))) + exp(data_j$Y[i]*eta_y[2]-log(1+exp(eta_y[2]))+eta_c)
     first.deriv[[i]] <- (exp(data_j$Y[i]*eta_y[1]-log(1+exp(eta_y[1])))*(data_j$Y[i]-Py[1])* deriv.etay.lambda[1] +
                            exp(data_j$Y[i]*eta_y[2]-log(1+exp(eta_y[2]))+eta_c)*(data_j$Y[i]-Py[2])*deriv.etay.lambda[2]) / D00
    }
    
    if(data_j$D[i]==1 & data_j$trt[i]==1){ #trt complier
     first.deriv[[i]] <- (data_j$Y[i]-Py[3]) * deriv.etay.lambda[3]
    }

  }#end of loop
  
  first_der_lj_lambda <- Reduce('+', first.deriv)

  return(first_der_lj_lambda) 
}






