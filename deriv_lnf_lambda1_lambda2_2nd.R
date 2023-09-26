#fun_2nd_deriv_lnf_lambda1_lambda2() gives the second derivative of lnf / lambda_k*lambda_l
#inputs:
#u: random effects of y model. usually the vector z of length r' after change of variable
#b: random effect of c model. usually the last element of vector z
#Alpha: alpha vector in Y model
#gamma: gamma vector in C model
#data_j: data frame of jth site
#ind.x0: column number of covariates controlled in Y model
#ind.x1: column number of covariates controlled in C model
#lambda: factor loading matrix
#M is the direvative of factor loading matrix wrt lambda_k and lambda_l

fun_2nd_deriv_lnf_lambda1_lambda2 <- function(u,b,Alpha,gamma,data_j,ind.x0,ind.x1,lambda,M){
  data_j <- data_j[!is.na(data_j$Y),]
  n.j <- nrow(data_j)
  
  if(anyNA(ind.x1)){
    Xij <- as.matrix(rep(1,n.j)) #intercept only for compliance model
  } else {
    Xij <- as.matrix(cbind(rep(1,n.j),data_j[,ind.x1])) #intercept and covariates
  }

  Bij <- diag(rep(1,3))

  deriv_2nd_lnf_lambda1_lambda2 <- list()
  
  for(i in 1:nrow(data_j)){
    if(anyNA(ind.x0)){
      Aij <- diag(rep(1,3))
    } else {
      Aij <- cbind(diag(rep(1,3)),matrix(unlist(replicate(3,as.matrix(data_j[,ind.x0])[i,])),nrow = 3,byrow = T))
    }
    
    eta_y <- Aij %*% Alpha + Bij %*% lambda %*% u
    Py <- 1/(1+exp(-eta_y)) #probability of y being 1. a vector of 3.
    eta_c <- as.numeric(t(as.matrix(gamma)) %*% as.matrix(Xij[i,]) + b)
    
    if(r_prime==1){
      deriv.etay.lambda1 <- t(Bij) %*% M[,1] %*% u
      deriv.etay.lambda2 <- t(Bij) %*% M[,2] %*% u
      
      if(data_j$D[i]==0 & data_j$trt[i]==1){ #never taker
        deriv_2nd_lnf_lambda1_lambda2[[i]] <- -Py[1]*(1-Py[1]) * deriv.etay.lambda1[1]*deriv.etay.lambda2[1]
      }
      
      if(data_j$D[i]==0 & data_j$trt[i]==0){ #control complier+never taker
        D00 <- exp(data_j$Y[i]*eta_y[1]-log(1+exp(eta_y[1]))) + exp(data_j$Y[i]*eta_y[2]-log(1+exp(eta_y[2]))+eta_c)
        
        D00_p_lambda_l <- exp(data_j$Y[i]*eta_y[1]-log(1+exp(eta_y[1])))*(data_j$Y[i]-Py[1]) * deriv.etay.lambda1[1] +
          exp(data_j$Y[i]*eta_y[2]-log(1+exp(eta_y[2]))+eta_c)*(data_j$Y[i]-Py[2]) * deriv.etay.lambda1[2]
        
        D00_p_lambda_k <- exp(data_j$Y[i]*eta_y[1]-log(1+exp(eta_y[1])))*(data_j$Y[i]-Py[1]) * deriv.etay.lambda2[1] +
          exp(data_j$Y[i]*eta_y[2]-log(1+exp(eta_y[2]))+eta_c)*(data_j$Y[i]-Py[2]) * deriv.etay.lambda2[2]
        
        D00_pp_lambda <- exp(data_j$Y[i]*eta_y[1]-log(1+exp(eta_y[1])))* deriv.etay.lambda1[1] * deriv.etay.lambda2[1] * ((data_j$Y[i]-Py[1])^2- Py[1]*(1-Py[1]))+
          exp(data_j$Y[i]*eta_y[2]-log(1+exp(eta_y[2]))+eta_c)* deriv.etay.lambda1[2] * deriv.etay.lambda2[2] * ((data_j$Y[i]-Py[2])^2- Py[2]*(1-Py[2]))
        
        deriv_2nd_lnf_lambda1_lambda2[[i]] <- ( D00_pp_lambda*D00 - D00_p_lambda_l %*% D00_p_lambda_k )/D00^2
      }
      
      if(data_j$D[i]==1 & data_j$trt[i]==1){ #trt complier
        deriv_2nd_lnf_lambda1_lambda2[[i]] <- -Py[3]*(1-Py[3]) * deriv.etay.lambda1[3] * deriv.etay.lambda2[3]
      }
    } else {
      deriv.etay.lambda1 <- t(Bij) %*% M[,1:2] %*% u
      deriv.etay.lambda2 <- t(Bij) %*% M[,3:4] %*% u
      
      if(data_j$D[i]==0 & data_j$trt[i]==1){ #never taker
        deriv_2nd_lnf_lambda1_lambda2[[i]] <- -Py[1]*(1-Py[1]) * deriv.etay.lambda1[1,] %*% deriv.etay.lambda2[1,]
      }
      
      if(data_j$D[i]==0 & data_j$trt[i]==0){ #control complier+never taker
        D00 <- exp(data_j$Y[i]*eta_y[1]-log(1+exp(eta_y[1]))) + exp(data_j$Y[i]*eta_y[2]-log(1+exp(eta_y[2]))+eta_c)
        
        D00_p_lambda_l <- exp(data_j$Y[i]*eta_y[1]-log(1+exp(eta_y[1])))*(data_j$Y[i]-Py[1]) * deriv.etay.lambda1[1,] +
          exp(data_j$Y[i]*eta_y[2]-log(1+exp(eta_y[2]))+eta_c)*(data_j$Y[i]-Py[2]) * deriv.etay.lambda1[2,]
        
        D00_p_lambda_k <- exp(data_j$Y[i]*eta_y[1]-log(1+exp(eta_y[1])))*(data_j$Y[i]-Py[1]) * deriv.etay.lambda2[1,] +
          exp(data_j$Y[i]*eta_y[2]-log(1+exp(eta_y[2]))+eta_c)*(data_j$Y[i]-Py[2]) * deriv.etay.lambda2[2,]
        
        D00_pp_lambda <- exp(data_j$Y[i]*eta_y[1]-log(1+exp(eta_y[1])))* deriv.etay.lambda1[1,] %*% deriv.etay.lambda2[1,] * ((data_j$Y[i]-Py[1])^2- Py[1]*(1-Py[1]))+
          exp(data_j$Y[i]*eta_y[2]-log(1+exp(eta_y[2]))+eta_c)* deriv.etay.lambda1[2,] %*% deriv.etay.lambda2[2,] * ((data_j$Y[i]-Py[2])^2- Py[2]*(1-Py[2]))
        
        deriv_2nd_lnf_lambda1_lambda2[[i]] <- ( D00_pp_lambda*D00 - D00_p_lambda_l %*% D00_p_lambda_k )/D00^2
      }
      
      if(data_j$D[i]==1 & data_j$trt[i]==1){ #trt complier
        deriv_2nd_lnf_lambda1_lambda2[[i]] <- -Py[3]*(1-Py[3]) * deriv.etay.lambda1[3,] %*% deriv.etay.lambda2[3,]
      }
    }
    
  }#end of loop
  
  deriv_2nd_lnf_lambda1_lambda2 <- Reduce('+', deriv_2nd_lnf_lambda1_lambda2)

  return(deriv_2nd_lnf_lambda1_lambda2)
}

##fun_2nd_deriv_L_lambda1_lambda2() gives the second derivative Lj/lambda_k*lambda_l
fun_2nd_deriv_L_lambda1_lambda2 <- function(u,b,Alpha,gamma,data_j,ind.x0,ind.x1,lambda,M){
  
  deriv_2nd <- fun_2nd_deriv_lnf_lambda1_lambda2(u,b,Alpha,gamma,data_j,ind.x0,ind.x1,lambda,M) + fun_der_lnf_lambda(u,b,Alpha,gamma,data_j,ind.x0,ind.x1,lambda,M[,1]) %*% t(fun_der_lnf_lambda(u,b,Alpha,gamma,data_j,ind.x0,ind.x1,lambda,M[,2]))

  return(deriv_2nd)
}

