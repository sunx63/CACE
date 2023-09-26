
#deriv_lnf_rc_1st() gives the first derivative wrt r
#inputs:
#u: random effects of y model. usually the vector z of length r' after change of variable
#b: random effect of c model. usually the last element of vector z
#Alpha: alpha vector in Y model
#gamma: gamma vector in C model
#data_j: data frame of jth site
#ind.x0: column number of covariates controlled in Y model
#ind.x1: column number of covariates controlled in C model
#lambda: factor loading matrix
deriv_lnf_rc_1st <- function(u,b,Alpha,gamma,data_j,ind.x0,ind.x1,lambda){
  data_j <- data_j[!is.na(data_j$Y) | !is.na(data_j$C),] #remove obs with both missing C and Y
  
  n.j <- nrow(data_j)

  if(anyNA(ind.x1)){
    Xij <- as.matrix(rep(1,n.j)) #intercept only for compliance model
  } else {
    Xij <- as.matrix(cbind(rep(1,n.j),data_j[,ind.x1])) #intercept and covariates
  }
  
  Bij <- diag(rep(1,3))
  Bij <- Bij %*% lambda
  
  first <- 0
  second <- 0
  third <- 0
  forth <- 0
  fifth <- 0
  
  for(i in 1:n.j){
    if(is.na(data_j$Y[i])){
      eta_c <- as.numeric(t(as.matrix(gamma)) %*% as.matrix(Xij[i,]) + b) #compliance model
      Pc <- 1/(1+exp(-eta_c))
      Pn <- 1-Pc
      if(data_j$trt[i]==1 & data_j$D[i]==0) forth <- forth - Pc*Xij[i,]
      if(data_j$trt[i]==1 & data_j$D[i]==1) fifth <- fifth + Pn*Xij[i,]
      
    } else {
      
    if(anyNA(ind.x0)){
      Aij <- diag(rep(1,3))
    } else {
      Aij <- cbind(diag(rep(1,3)),matrix(unlist(replicate(3,as.matrix(data_j[,ind.x0])[i,])),nrow = 3,byrow = T))
    }
    
    eta_y <- Aij %*% Alpha + Bij %*% u #y model
    eta_c <- as.numeric(t(as.matrix(gamma)) %*% as.matrix(Xij[i,]) + b)
    
    if(data_j$trt[i]==0 & data_j$D[i]==0){
      first <- first + exp(data_j$Y[i]*eta_y[2]-log(1+exp(eta_y[2]))+eta_c) * Xij[i,] /
        (exp(data_j$Y[i]*eta_y[1]-log(1+exp(eta_y[1])))+exp(data_j$Y[i]*eta_y[2]-log(1+exp(eta_y[2]))+eta_c))
    }
    if(data_j$trt[i]==1 & data_j$D[i]==1){
      second <- second + as.matrix(Xij[i,])
    } #end of if
    third <- third + as.matrix(exp(eta_c) * Xij[i,] / (1+exp(eta_c)))
    }#end of else
  } #end of loop
  
  return(first + second - third + forth + fifth)
}#end of function




