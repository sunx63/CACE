
#f_y_incomplete() calculates the joint pmf of Y and C with C partially observed.
#inputs:
#u: random effects of y model. usually the vector z of length r' after change of variable
#b: random effect of c model. usually the last element of vector z
#Alpha: alpha vector in Y model
#gamma: gamma vector in C model
#data_j: data frame of jth site
#ind.x0: column number of covariates controlled in Y model
#ind.x1: column number of covariates controlled in C model
#lambda: factor loading matrix

f_y_incomplete <- function(u,b,Alpha,gamma,data_j,ind.x0,ind.x1,lambda){
   data_j <- data_j[!is.na(data_j$Y) | !is.na(data_j$C),] #remove obs with both missing C and Y

   n.j <- nrow(data_j)
  
   if(anyNA(ind.x1)){
     Xij <- as.matrix(rep(1,n.j)) #intercept only for compliance model
   } else {
     Xij <- as.matrix(cbind(rep(1,n.j), data_j[,ind.x1])) #intercept and covariates
   }
   
   Bij <- diag(rep(1,3))
   Bij <- Bij %*% lambda

   data_j$f_ij <- 0

   for(i in 1:n.j){
     if(is.na(data_j$Y[i])){
       eta_c <- as.numeric(t(as.matrix(gamma)) %*% as.matrix(Xij[i,]) + b) #compliance model. eta_c
       Pc <- 1/(1+exp(-eta_c))
       Pn <- 1-Pc
       if(data_j$trt[i]==1 & data_j$D[i]==0) data_j$f_ij[i] <- Pn
       if(data_j$trt[i]==1 & data_j$D[i]==1) data_j$f_ij[i] <- Pc
     } else {
       if(anyNA(ind.x0)){
         Aij <- diag(rep(1,3))
       } else {
         Aij <- cbind(diag(rep(1,3)),matrix(unlist(replicate(3,as.matrix(data_j[,ind.x0])[i,])),nrow = 3,byrow = T))
       }
       
       eta_y <- Aij %*% as.matrix(Alpha) + Bij %*% u
       eta_c <- as.numeric(t(as.matrix(gamma)) %*% as.matrix(Xij[i,]) + b) 
       
       Pc <- 1/(1+exp(-eta_c))
       Pn <- 1-Pc
       
       if(data_j$trt[i]==1 & data_j$D[i]==0){ #never taker
         data_j$f_ij[i] <- exp(data_j$Y[i]*eta_y[1]-log(1+exp(eta_y[1])))*Pn 
       }
       if(data_j$trt[i]==0 & data_j$D[i]==0){ #control complier+never taker
         data_j$f_ij[i] <- exp(data_j$Y[i]*eta_y[1]-log(1+exp(eta_y[1])))*Pn + exp(data_j$Y[i]*eta_y[2]-log(1+exp(eta_y[2])))*Pc
       }
       if(data_j$trt[i]==1 & data_j$D[i]==1){ #trt complier 
         data_j$f_ij[i] <- exp(data_j$Y[i]*eta_y[3]-log(1+exp(eta_y[3])))*Pc 
       }
     }#end of ifelse
   }#end of loop
   if(0 %in% data_j$f_ij){
     f_ij <- data_j$f_ij[-which(data_j$f_ij %in% 0)]
   } else {
     f_ij <- data_j$f_ij
   }
    
    f_y_j <- prod(f_ij)
    
   return(f_y_j)
 }


