
#Main functions for estimating CACE with parallel computing
#inputs
#r: full dimension of random effects of Y model, numeric. r=3 for one-sided noncompliance
#k: full dimension of random effects of C model, numeric. k=1 for one-sided noncompliance
#Q: number of abscissas for AGHQ, numeric
#J: number of clusters/sites, numeric
#x0: covariates controlled in Y model, character
#x1: covariates controlled in C model, character
#init: initial values for parameters, numeric
#data: dataframe containing outcome, treatment assignment and receipt, cluster/site level variable 
#col.clinic: column number of cluster/site level variable in the dataframe, numeric
#col.trt: column number of treatment assignment in the dataframe, numeric
#col.D: column number of treatment receipt in the dataframe, numeric
#col.Y: column number of outcome variable in the dataframe, numeric
#niter: maximum number of iteration, numeric. suggest 500
#tol: convergence criteria, numeric. suggest 10^-4
#Share: sharing of random effect, numeric. 
  #When r=1, Share=1 means u'=un, Share=2 means u'=uc0 and Share=3 means u'=uc1
  #When r=2, Share=1 means never takers and treatment compliers sharing one random effect, Share=2 means never takers and control compliers sharing one random effect, Share=3 means control compliers and trt compliers sharing one random effect
#r_prime: reduced number of random effects, numeric

cace_1side_parallel_lambda <- function(r=3,k=1,Q,J,data,x0,x1,init,col.clinic,col.trt,col.D,col.Y,niter,tol,Share,r_prime){ 
  start_time <- Sys.time()
  
  #first derivative lnf/a * f(y)
  lnf_a_fy_1 <- function(u,b,Alpha,gamma,data_j,ind.x0,ind.x1,lambda){
    fun_der_lnf_alpha(u,b,Alpha,gamma,data_j,ind.x0,ind.x1,lambda) * f_y_incomplete(u,b,Alpha,gamma,data_j,ind.x0,ind.x1,lambda)
  }
  
  #second derivative lnf/a * f(y)
  lnf_a_fy_2 <- function(u,b,Alpha,gamma,data_j,ind.x0,ind.x1,lambda){
    fun_2nd_deriv_L_a(u,b,Alpha,gamma,data_j,ind.x0,ind.x1,lambda) * f_y_incomplete(u,b,Alpha,gamma,data_j,ind.x0,ind.x1,lambda)
  } 
  
  #first derivative lnf/tau * f(y)
  lnf_t_fy_1 <- function(u,b,Alpha,gamma,data_j,ind.x0,ind.x1,lambda,tau){
    fun_der_1st_ln_phitau(u,tau) * f_y_incomplete(u,b,Alpha,gamma,data_j,ind.x0,ind.x1,lambda)
  }
  
  #second derivative lnf/tau * f(y)
  lnf_t_fy_2 <- function(u,b,Alpha,gamma,data_j,ind.x0,ind.x1,lambda,tau){
    fun_2nd_L_tau(u,tau) * f_y_incomplete(u,b,Alpha,gamma,data_j,ind.x0,ind.x1,lambda)
  }
  
  #first derivative lnf/rc * f(y)
  lnf_rc_fy_1 <- function(u,b,Alpha,gamma,data_j,ind.x0,ind.x1,lambda){
    deriv_lnf_rc_1st(u,b,Alpha,gamma,data_j,ind.x0,ind.x1,lambda) * f_y_incomplete(u,b,Alpha,gamma,data_j,ind.x0,ind.x1,lambda)
  }
  
  #seconc derivative lnf/rc * f(y)
  lnf_rc_fy_2 <- function(u,b,Alpha,gamma,data_j,ind.x0,ind.x1,lambda){
    fun_2nd_deriv_L_rc(u,b,Alpha,gamma,data_j,ind.x0,ind.x1,lambda) * f_y_incomplete(u,b,Alpha,gamma,data_j,ind.x0,ind.x1,lambda)
  } 
  
  #second derivative lnf/tau*a * f(y)
  lnf_t_a_fy_2 <- function(u,b,Alpha,gamma,data_j,ind.x0,ind.x1,lambda,tau){
    f_y_incomplete(u,b,Alpha,gamma,data_j,ind.x0,ind.x1,lambda) * fun_der_1st_ln_phitau(u,tau) %*% fun_der_lnf_alpha(u,b,Alpha,gamma,data_j,ind.x0,ind.x1,lambda) 
  }
  
  #first derivative lnf/delta * f(y)
  lnf_delta_fy_1 <- function(u,b,Alpha,gamma,data_j,ind.x0,ind.x1,lambda,delta){
    deriv_lnf_delta_1st(b,delta) * f_y_incomplete(u,b,Alpha,gamma,data_j,ind.x0,ind.x1,lambda)
  }
  
  #second derivative lnf/delta * f(y)
  lnf_delta_fy_2 <- function(u,b,Alpha,gamma,data_j,ind.x0,ind.x1,lambda,delta){
    fun_2nd_L_delta(b,delta) * f_y_incomplete(u,b,Alpha,gamma,data_j,ind.x0,ind.x1,lambda)
  }
  
  #second derivative lnf/delta*a * f(y)
  lnf_delta_a_fy_2 <- function(u,b,Alpha,gamma,data_j,ind.x0,ind.x1,lambda,delta){
    f_y_incomplete(u,b,Alpha,gamma,data_j,ind.x0,ind.x1,lambda) * deriv_lnf_delta_1st(b,delta) %*% fun_der_lnf_alpha(u,b,Alpha,gamma,data_j,ind.x0,ind.x1,lambda) 
  }
  
  #second derivative lnf/tau*rc * f(y)
  lnf_t_rc_fy_2 <- function(u,b,Alpha,gamma,data_j,ind.x0,ind.x1,lambda,tau){
    f_y_incomplete(u,b,Alpha,gamma,data_j,ind.x0,ind.x1,lambda) * fun_der_1st_ln_phitau(u,tau) %*% t(deriv_lnf_rc_1st(u,b,Alpha,gamma,data_j,ind.x0,ind.x1,lambda))
  }
  
  #second derivative lnf/delta*rc * f(y)
  lnf_delta_rc_fy_2 <- function(u,b,Alpha,gamma,data_j,ind.x0,ind.x1,lambda,delta){
    f_y_incomplete(u,b,Alpha,gamma,data_j,ind.x0,ind.x1,lambda) * deriv_lnf_delta_1st(b,delta) %*% t(deriv_lnf_rc_1st(u,b,Alpha,gamma,data_j,ind.x0,ind.x1,lambda))
  }
  
  #second derivative lnf/tau*delta * f(y)
  lnf_delta_t_fy_2 <- function(u,b,Alpha,gamma,data_j,ind.x0,ind.x1,lambda,tau,delta){
    deriv_lnf_delta_1st(b,delta) %*% t(fun_der_1st_ln_phitau(u,tau)) * f_y_incomplete(u,b,Alpha,gamma,data_j,ind.x0,ind.x1,lambda)
  } 
  
  #second derivative lnf/rc*a * f(y)
  lnf_rc_a_fy_2 <- function(u,b,Alpha,gamma,data_j,ind.x0,ind.x1,lambda){
    f_y_incomplete(u,b,Alpha,gamma,data_j,ind.x0,ind.x1,lambda)*fun_2nd_deriv_L_alpha_rc(u,b,Alpha,gamma,data_j,ind.x0,ind.x1,lambda)
  }
  
  #f(y)*z E(u,b)
  f_y_z <- function(u,b,Alpha,gamma,data_j,ind.x0,ind.x1,lambda){
    f_y_incomplete(u,b,Alpha,gamma,data_j,ind.x0,ind.x1,lambda) * c(u,b)
  }
  
  #f(y)*z^2 E((u,b)^2)
  f_y_z2 <- function(u,b,Alpha,gamma,data_j,ind.x0,ind.x1,lambda){
    f_y_incomplete(u,b,Alpha,gamma,data_j,ind.x0,ind.x1,lambda) * c(u,b) %*% t(c(u,b))
  }
  
  #first derivative lnf/lambda * f(y)
  lnf_lambda_fy_1 <- function(u,b,Alpha,gamma,data_j,ind.x0,ind.x1,lambda,M){
    fun_der_lnf_lambda(u,b,Alpha,gamma,data_j,ind.x0,ind.x1,lambda,M) * f_y_incomplete(u,b,Alpha,gamma,data_j,ind.x0,ind.x1,lambda)
  }
  
  #second derivative lnf/lambda * f(y)
  lnf_lambda_fy_2 <- function(u,b,Alpha,gamma,data_j,ind.x0,ind.x1,lambda,M){
    fun_2nd_deriv_L_lambda(u,b,Alpha,gamma,data_j,ind.x0,ind.x1,lambda,M) * f_y_incomplete(u,b,Alpha,gamma,data_j,ind.x0,ind.x1,lambda)
  }
  
  #second derivative lnf/a*lambda * f(y)
  lnf_a_lambda_fy_2 <- function(u,b,Alpha,gamma,data_j,ind.x0,ind.x1,lambda,M){
    fun_2nd_deriv_L_lambda_a(u,b,Alpha,gamma,data_j,ind.x0,ind.x1,lambda,M) * f_y_incomplete(u,b,Alpha,gamma,data_j,ind.x0,ind.x1,lambda)
  }
  
  #second derivative lnf/rc*lambda * f(y)
  lnf_rc_lambda_fy_2 <- function(u,b,Alpha,gamma,data_j,ind.x0,ind.x1,lambda,M){
    fun_2nd_deriv_L_lambda_rc(u,b,Alpha,gamma,data_j,ind.x0,ind.x1,lambda,M) * f_y_incomplete(u,b,Alpha,gamma,data_j,ind.x0,ind.x1,lambda)
  }
  
  #second derivative lnf/lambda1*lambda2 * f(y)
  lnf_lambda1_lambda2_fy_2 <- function(u,b,Alpha,gamma,data_j,ind.x0,ind.x1,lambda,M){
    fun_2nd_deriv_L_lambda1_lambda2(u,b,Alpha,gamma,data_j,ind.x0,ind.x1,lambda,M) * f_y_incomplete(u,b,Alpha,gamma,data_j,ind.x0,ind.x1,lambda)
  }
  
  
  lnf_t_lambda_fy_2 <- function(u,b,Alpha,gamma,data_j,ind.x0,ind.x1,lambda,tau,M){
    f_y_incomplete(u,b,Alpha,gamma,data_j,ind.x0,ind.x1,lambda) * fun_der_1st_ln_phitau(u,tau) %*% fun_der_lnf_lambda(u,b,Alpha,gamma,data_j,ind.x0,ind.x1,lambda,M) 
  }
  
  lnf_delta_lambda_fy_2 <- function(u,b,Alpha,gamma,data_j,ind.x0,ind.x1,lambda,delta,M){
    f_y_incomplete(u,b,Alpha,gamma,data_j,ind.x0,ind.x1,lambda) * deriv_lnf_delta_1st(b,delta) %*% fun_der_lnf_lambda(u,b,Alpha,gamma,data_j,ind.x0,ind.x1,lambda,M) 
  }
  
  colnames(data)[c(col.clinic,col.trt, col.D, col.Y)] <- c("clinic", "trt", "D", "Y") 
  
  if(anyNA(x0)){
    ind.x0 <- NA
    l.alpha <- 3
  } else {
    ind.x0 <- which(colnames(data) %in% c(x0))
    l.alpha <- 3 + length(x0)
  }
  
  if(anyNA(x1)){
    ind.x1 <- NA
    l.gamma <- 1
  } else {
    ind.x1 <- which(colnames(data) %in% c(x1)) #x1 should be input as a vector of character with "". covariate order need to comply with the order of collumns in dataset
    l.gamma <- 1+length(x1)
  }
  
  for(i in 1:nrow(data)){
    if(data$trt[i]==1 & data$D[i]==1){
      data$C[i] <- 1
    } 
    if(data$trt[i]==1 & data$D[i]==0){
      data$C[i] <- 0
    }
    if(data$trt[i]==0 & data$D[i]==0){
      data$C[i] <- NA
    }
  }
  
  nlambda <- r-r_prime
  
  if(nlambda==1){
    lambda1 <- init[l.alpha+l.gamma+1]
  }
  if(nlambda==2){
    lambda1 <- init[l.alpha+l.gamma+1]
    lambda2 <- init[l.alpha+l.gamma+2]
  }
 
  if(!missing(r_prime)){
    if(r_prime==1){ 
      if(Share==1) lambda <- as.matrix(c(1,lambda1,lambda2)) #share never taker's random effect. u'=un
      if(Share==2) lambda <- as.matrix(c(lambda1,1,lambda2)) #share control complier's random effect. u'=uc0
      if(Share==3) lambda <- as.matrix(c(lambda1,lambda2,1)) ##share trt complier's random effect. u'=uc1
    }
    if(r_prime==2){ 
      if(Share==1) lambda <- matrix(c(1,0,lambda1,0,1,0),nrow = 3,ncol = r_prime,byrow = F) #never takers and treatment compliers sharing one random effect
      if(Share==2) lambda <- matrix(c(1,lambda1,0,0,0,1),nrow = 3,ncol = r_prime,byrow = F) #never takers and control compliers sharing one random effect
      if(Share==3) lambda <- matrix(c(1,0,0,0,1,lambda1),nrow = 3,ncol = r_prime,byrow = F) #control compliers and trt compliers sharing one random effect
    }
  }

  theta <- matrix(NA,niter,(r_prime*(r_prime+1)/2+l.alpha+l.gamma+k*(k+1)/2+nlambda)) 
  #  'alpha_n','alpha_c0','alpha_c1',gamma,lambda1,lambda2,'tou11','tou12','tou13','tou22','tou23','tou33',delta)
  theta[1,] <- init 
  
  if(r_prime==1){
    Vu <- theta[1,(l.alpha+l.gamma+nlambda+1):(r_prime*(r_prime+1)/2+l.alpha+l.gamma+nlambda)] 
  } else {
    Vu <- VechToCovM(theta[1,(l.alpha+l.gamma+nlambda+1):(r_prime*(r_prime+1)/2+l.alpha+l.gamma+nlambda)], r_prime)
  }
  
  if(k==1){
    Vb <- init[((length(init)-k*(k+1)/2+1):length(init))]
  } else {
    Vb <- VechToCovM(init[((length(init)-k*(k+1)/2+1):length(init))], k)
  }
  

  Vub <- bdiag(Vu, Vb) %>% as.matrix
  Vu.b <- do.call("rbind",replicate(J, MtoVech(Vub,(r_prime+k)), simplify = FALSE))
  
  E_u <- do.call("rbind",replicate(J, rep(0,(r_prime+k)), simplify = FALSE))
  
  E_u2 <- NA
  
  X <- gaussHermiteData(Q)
  
  #aghq matrix for u
  A <- matrix(NA, Q^(r_prime+k),(r_prime+k))
  
  for(i in 1:(r_prime+k)){
    A[,i] <- rep(X$x,each=Q^((r_prime+k)-i))
  }
  
  W <- matrix(NA, Q^(r_prime+k),(r_prime+k))
  for(i in 1:(r_prime+k)){
    W[,i] <- rep(X$w,each=Q^((r_prime+k)-i))
  }
  
  Ws <- W*exp(A^2)
  Ws <- apply(Ws,1,prod)
  
  L2.norm <- c() #l2 norm between theta^k and theta^(k+1)
  ll <- c() #likelihood
  
  if(r_prime==1){
    M1 <- c(0,1,0)
    M2 <- c(0,0,1)
    M12 <- cbind(M1,M2)
  }
  if(r_prime==2){
    if(Share==1) M1 <- matrix(c(0,0,1,rep(0,3)),nrow=3,ncol=2,byrow = F)
    if(Share==2) M1 <- matrix(c(0,1,rep(0,4)),nrow=3,ncol=2,byrow = F)
    if(Share==3) M1 <- matrix(c(rep(0,5),1),nrow=3,ncol=2,byrow = F)
  }
  
  for(iter in 1:niter){
    Alpha <- theta[iter,1:l.alpha]
    gamma <- theta[iter,(l.alpha+1):(l.alpha+l.gamma)]
    
    if(r_prime==1){
      tau <- theta[iter,(length(c(Alpha,gamma))+nlambda+1):(r_prime*(r_prime+1)/2+l.alpha+l.gamma+nlambda)]
    } else {
      tau <- VechToCovM(theta[iter,(length(c(Alpha,gamma))+nlambda+1):(r_prime*(r_prime+1)/2+l.alpha+l.gamma+nlambda)],r_prime)
    }
    
    if(k==1){
      delta <- theta[iter,((length(init)-k*(k+1)/2+1):length(init))]
    } else {
      delta <- VechToCovM(theta[iter,((length(init)-k*(k+1)/2+1):length(init))],k)
    }
    
    tau.delta <- bdiag(tau, delta) %>% as.matrix
    
    print(paste0("iteration = ", iter))
    
    # Integration by AGHQ starts here
    list1 <- foreach(j = 1:J,.verbose=T,.errorhandling='stop',.export=ls(.GlobalEnv),.packages=c('mvtnorm')) %dopar% { 
      #filter jth clinic data
      data_j <- data[data$clinic==j,]
      
      #change of variable
      V.z <- VechToCovM(Vu.b[j,],(k+r_prime))
      
      Lz <- t(chol(2*V.z))
      z <- matrix(NA,Q^(r_prime+k),(r_prime+k))
      
      Wphi <- numeric(length=nrow(z))
      for(i in 1:nrow(z)){
        z[i,] <- Lz %*% A[i,]+E_u[j,]
        Wphi[i] <- Ws[i]*dmvnorm(z[i,], mean = rep(0, (r_prime+k)), sigma = tau.delta, log = FALSE)
      }
      
      #AGHQ     
      Lj <- det(Lz)*AGHQ(Q,Wphi,f_y_incomplete,z,r_prime,k,Alpha,gamma,data_j,ind.x0,ind.x1,lambda)
      if(is.nan(Lj) | Lj==0) stop("lj is not positive")
      
      E_u_j <- (det(Lz)/Lj) * AGHQ(Q,Wphi,f_y_z,z,r_prime,k,Alpha,gamma,data_j,ind.x0,ind.x1,lambda)
      E_u2 <- (det(Lz)/Lj) * AGHQ(Q,Wphi,f_y_z2,z,r_prime,k,Alpha,gamma,data_j,ind.x0,ind.x1,lambda)
      
      V <- E_u2 - E_u_j %*% t(E_u_j)

      Vub_j <- MtoVech(V,(r_prime+k))
      
      S_L_a <- det(Lz) * AGHQ(Q,Wphi,lnf_a_fy_1,z,r_prime,k,Alpha,gamma,data_j,ind.x0,ind.x1,lambda) 
      H_L_a <- det(Lz) * AGHQ(Q,Wphi,lnf_a_fy_2,z,r_prime,k,Alpha,gamma,data_j,ind.x0,ind.x1,lambda) 
      S_L_rc <- det(Lz) * AGHQ(Q,Wphi,lnf_rc_fy_1,z,r_prime,k,Alpha,gamma,data_j,ind.x0,ind.x1,lambda)
      H_L_rc <- det(Lz) * AGHQ(Q,Wphi,lnf_rc_fy_2,z,r_prime,k,Alpha,gamma,data_j,ind.x0,ind.x1,lambda)
      H_L_a_rc <- det(Lz) * AGHQ(Q,Wphi,lnf_rc_a_fy_2,z,r_prime,k,Alpha,gamma,data_j,ind.x0,ind.x1,lambda)

      
      if(nlambda==1){
        S_L_lambda1 <- det(Lz) * AGHQ(Q,Wphi,lnf_lambda_fy_1,z,r_prime,k,Alpha,gamma,data_j,ind.x0,ind.x1,lambda,M1)
        H_L_lambda1 <- det(Lz) * AGHQ(Q,Wphi,lnf_lambda_fy_2,z,r_prime,k,Alpha,gamma,data_j,ind.x0,ind.x1,lambda,M1)
        H_L_a_lambda1 <- det(Lz) * AGHQ(Q,Wphi,lnf_a_lambda_fy_2,z,r_prime,k,Alpha,gamma,data_j,ind.x0,ind.x1,lambda,M1)
        H_L_rc_lambda1 <- det(Lz) * AGHQ(Q,Wphi,lnf_rc_lambda_fy_2,z,r_prime,k,Alpha,gamma,data_j,ind.x0,ind.x1,lambda,M1)
      }
      
      if(nlambda==2){
        S_L_lambda1 <- det(Lz) * AGHQ(Q,Wphi,lnf_lambda_fy_1,z,r_prime,k,Alpha,gamma,data_j,ind.x0,ind.x1,lambda,M1)
        S_L_lambda2 <- det(Lz) * AGHQ(Q,Wphi,lnf_lambda_fy_1,z,r_prime,k,Alpha,gamma,data_j,ind.x0,ind.x1,lambda,M2)
        H_L_lambda1 <- det(Lz) * AGHQ(Q,Wphi,lnf_lambda_fy_2,z,r_prime,k,Alpha,gamma,data_j,ind.x0,ind.x1,lambda,M1)
        H_L_lambda2 <- det(Lz) * AGHQ(Q,Wphi,lnf_lambda_fy_2,z,r_prime,k,Alpha,gamma,data_j,ind.x0,ind.x1,lambda,M2)
        H_L_a_lambda1 <- det(Lz) * AGHQ(Q,Wphi,lnf_a_lambda_fy_2,z,r_prime,k,Alpha,gamma,data_j,ind.x0,ind.x1,lambda,M1)
        H_L_a_lambda2 <- det(Lz) * AGHQ(Q,Wphi,lnf_a_lambda_fy_2,z,r_prime,k,Alpha,gamma,data_j,ind.x0,ind.x1,lambda,M2)
        H_L_rc_lambda1 <- det(Lz) * AGHQ(Q,Wphi,lnf_rc_lambda_fy_2,z,r_prime,k,Alpha,gamma,data_j,ind.x0,ind.x1,lambda,M1)
        H_L_rc_lambda2 <- det(Lz) * AGHQ(Q,Wphi,lnf_rc_lambda_fy_2,z,r_prime,k,Alpha,gamma,data_j,ind.x0,ind.x1,lambda,M2)
        H_L_lambda1_lambda2 <- det(Lz) * AGHQ(Q,Wphi,lnf_lambda1_lambda2_fy_2,z,r_prime,k,Alpha,gamma,data_j,ind.x0,ind.x1,lambda,M12)
      }
 
      logLj <- log(Lj)
    
      Sa_j <- S_L_a/Lj
      Src_j <- S_L_rc/Lj
      
      Ha_j <- (H_L_a/Lj - S_L_a %*% t(S_L_a)/(Lj^2)) 
      Hrc_j <-  (H_L_rc/Lj - S_L_rc %*% t(S_L_rc)/(Lj^2)) 
      Harc_j <- (H_L_a_rc/Lj - S_L_rc %*% t(S_L_a)/(Lj^2)) 
      
      Slambda1_j <- S_L_lambda1/Lj
      Hlambda1_j <- (H_L_lambda1/Lj - S_L_lambda1 %*% t(S_L_lambda1)/(Lj^2))
      Halambda1_j <- (H_L_a_lambda1/Lj - S_L_lambda1 %*% t(S_L_a)/(Lj^2))
      Hrclambda1_j <- (H_L_rc_lambda1/Lj - S_L_lambda1 %*% t(S_L_rc)/(Lj^2))
      
      Slambda2_j <-0
      Hlambda2_j <-0
      Halambda2_j <-0
      Hrclambda2_j <-0
      Hlambda1_2_j <-0
      
      if(nlambda==2){
        Slambda2_j <- S_L_lambda2/Lj
        Hlambda2_j <- (H_L_lambda2/Lj - S_L_lambda2 %*% t(S_L_lambda2)/(Lj^2))
        Halambda2_j <- (H_L_a_lambda2/Lj - S_L_lambda2 %*% t(S_L_a)/(Lj^2))
        Hrclambda2_j <- (H_L_rc_lambda2/Lj - S_L_lambda2 %*% t(S_L_rc)/(Lj^2))
        Hlambda1_2_j <- (H_L_lambda1_lambda2/Lj - S_L_lambda2 %*% t(S_L_lambda1)/(Lj^2))
       
      }

      list2 <- list(E_u_j=E_u_j,E_u2=E_u2,Vub_j=Vub_j, logLj=logLj, Sa_j=Sa_j, Src_j=Src_j, Ha_j=Ha_j, Hrc_j=Hrc_j, Harc_j=Harc_j,
                    Slambda1_j=Slambda1_j, Hlambda1_j=Hlambda1_j, Halambda1_j=Halambda1_j, Hrclambda1_j=Hrclambda1_j,
                    Slambda2_j=Slambda2_j, Hlambda2_j=Hlambda2_j, Halambda2_j=Halambda2_j, Hrclambda2_j=Hrclambda2_j, Hlambda1_2_j=Hlambda1_2_j)
    }#end of loop for clinic
    
    Sa <- 0
    Src <- 0  
    
    Ha <- 0
    Harc <- 0
    Hrc <- 0
    
    Slambda1 <- 0
    Slambda2 <- 0
    Hlambda1 <- 0
    Hlambda2 <- 0
    Halambda1 <- 0
    Halambda2 <- 0
    Hrclambda1 <- Hrclambda2 <- 0
    Hlambda1_2 <- 0
    
    logL <- 0
    tau.del.hat <-0
 
    for(j in 1:J) {
      tau.del.hat <- tau.del.hat+list1[[j]]$E_u2
      logL <- logL+list1[[j]]$logLj
      Sa <- Sa+list1[[j]]$Sa_j
      Src <- Src+list1[[j]]$Src_j
      Ha <- Ha+list1[[j]]$Ha_j
      Hrc <- Hrc+list1[[j]]$Hrc_j
      Harc <- Harc+list1[[j]]$Harc_j
      Slambda1 <- Slambda1+list1[[j]]$Slambda1_j
      Slambda2 <- Slambda2+list1[[j]]$Slambda2_j
      Hlambda1 <- Hlambda1+list1[[j]]$Hlambda1_j
      Hlambda2 <- Hlambda2+list1[[j]]$Hlambda2_j
      Halambda1 <- Halambda1+list1[[j]]$Halambda1_j
      Halambda2 <- Halambda2+list1[[j]]$Halambda2_j
      Hrclambda1 <- Hrclambda1+list1[[j]]$Hrclambda1_j
      Hrclambda2 <- Hrclambda2+list1[[j]]$Hrclambda2_j
      Hlambda1_2 <- Hlambda1_2+list1[[j]]$Hlambda1_2_j
    }

    E_u <- list1[[1]]$E_u_j
    Vu.b <- list1[[1]]$Vub_j
    
    for(j in 2:J){
      E_u <- rbind(E_u, list1[[j]]$E_u_j)
      Vu.b <- rbind(Vu.b, list1[[j]]$Vub_j)
    }
    
    tau.d.h <- tau.del.hat/J
    tau.h <- tau.d.h[1:r_prime,1:r_prime]
    d.h <- tau.d.h[(nrow(tau.d.h)-k+1):nrow(tau.d.h),(nrow(tau.d.h)-k+1):nrow(tau.d.h)]
    #log likelihood of J clinics
    print(paste0("loglike = ", logL))
    
    if(nlambda==1){
      S <- c(Sa,Src,Slambda1)
   
      H1 <- cbind(Ha,t(Harc),t(Halambda1))
      H2 <- cbind(Harc,Hrc,t(Hrclambda1))
      H3 <- cbind(Halambda1,Hrclambda1,Hlambda1)
      H <- rbind(H1,H2,H3)
    }
    
    if(nlambda==2){
      S <- c(Sa,Src,Slambda1,Slambda2)
      
      H1 <- cbind(Ha,t(Harc),t(Halambda1),t(Halambda2))
      H2 <- cbind(Harc,Hrc,t(Hrclambda1),t(Hrclambda2))
      H3 <- cbind(Halambda1,Hrclambda1,Hlambda1,Hlambda1_2)
      H4 <- cbind(Halambda2,Hrclambda2,Hlambda1_2,Hlambda2)
      H <- rbind(H1,H2,H3,H4)
    }
    
    if(!is.symmetric.matrix(H)) {
      print(paste0("iter= ",iter))
      print(H)
      print("H asymmetric")
    }
    
    Hinv <- solve(H)
    
    theta[(iter+1),1:(l.alpha+l.gamma+nlambda)] <- theta[iter,1:(l.alpha+l.gamma+nlambda)] - as.vector(Hinv %*% S) 
    lambda1 <- theta[(iter+1),(l.alpha+l.gamma+1)]
    if(nlambda==2){
      lambda2 <- theta[(iter+1),(l.alpha+l.gamma+2)]
      if(Share==1) lambda <- c(1,lambda1,lambda2)
    }
    if(nlambda==1){
      if(Share==1)  lambda <- matrix(c(1,0,lambda1,0,1,0),nrow = 3,ncol = r_prime,byrow = F)
      if(Share==2)  lambda <- matrix(c(1,lambda1,rep(0,3),1),nrow = 3,ncol = r_prime, byrow = F)
      if(Share==3)  lambda <- matrix(c(1,rep(0,3),1,lambda1),nrow = 3, ncol = r_prime, byrow = F)
    }

    if(r_prime==1){
      theta[(iter+1),(l.alpha+l.gamma+1+nlambda):(r_prime*(r_prime+1)/2+l.alpha+l.gamma+nlambda)] <- tau.h
    } else {
      theta[(iter+1),(l.alpha+l.gamma+1+nlambda):(r_prime*(r_prime+1)/2+l.alpha+l.gamma+nlambda)] <- MtoVech(tau.h,r_prime)
    }
    if(k==1){
      theta[(iter+1),(l.alpha+l.gamma+r_prime*(r_prime+1)/2+nlambda+1):(r_prime*(r_prime+1)/2+l.alpha+l.gamma+k*(k+1)/2+nlambda)] <- d.h
    } else {
      theta[(iter+1),(l.alpha+l.gamma+r_prime*(r_prime+1)/2+nlambda+1):(r_prime*(r_prime+1)/2+l.alpha+l.gamma+k*(k+1)/2+nlambda)] <- MtoVech(d.h,k)
    }

    print(paste0("theta= ",theta[(iter+1),]))
    
    if(iter>1){
      L2.norm[iter] <- sqrt(sum((theta[(iter),]-theta[(iter-1),])^2))
      print(paste0("L2.norm= ",L2.norm[iter]))
      if(L2.norm[iter] < tol) {break}
    }
    
    ll[iter] <- logL
    
  } #end of iter loop
  
  para0 <- theta[(iter+1),]
  
  E_u2 <- NA
  
  Alpha <- para0[1:l.alpha]
  gamma <- para0[(l.alpha+1):(l.alpha+l.gamma)]
  
  if(r_prime==1){
    tau <- para0[(l.alpha+l.gamma+nlambda+1):(r_prime*(r_prime+1)/2+l.alpha+l.gamma+nlambda)]
  } else {
    tau <- VechToCovM(para0[(l.alpha+l.gamma+nlambda+1):(r_prime*(r_prime+1)/2+l.alpha+l.gamma+nlambda)], r_prime)
  }

  if(k==1){
    delta <- para0[((length(para0)-k*(k+1)/2+1):length(para0))]
  } else {
    delta <- VechToCovM(para0[((length(para0)-k*(k+1)/2+1):length(para0))], k)
  }
  
  tau.delta <- bdiag(tau, delta) %>% as.matrix


  #for computing H matrix only
  list1 <- foreach(j = 1:J,.verbose=T,.errorhandling='stop',.export=ls(.GlobalEnv),.packages=c('mvtnorm')) %dopar% { 
    #filter jth clinic data
    data_j <- data[data$clinic==j,]
    
    #change of variable
    V.z <- VechToCovM(Vu.b[j,],(k+r_prime))
    
    Lz <- t(chol(2*V.z))
    z <- matrix(NA,Q^(r_prime+k),(r_prime+k))
    
    Wphi <- numeric(length=nrow(z))
    for(i in 1:nrow(z)){
      z[i,] <- Lz %*% A[i,]+E_u[j,]
      Wphi[i] <- Ws[i]*dmvnorm(z[i,], mean = rep(0, (r_prime+k)), sigma = tau.delta, log = FALSE)
    }
    
    #AGHQ     
    Lj <- det(Lz)*AGHQ(Q,Wphi,f_y_incomplete,z,r_prime,k,Alpha,gamma,data_j,ind.x0,ind.x1,lambda)
    if(is.nan(Lj) | Lj==0) stop("lj is not positive")
    
    E_u_j <- (det(Lz)/Lj) * AGHQ(Q,Wphi,f_y_z,z,r_prime,k,Alpha,gamma,data_j,ind.x0,ind.x1,lambda)
    E_u2 <- (det(Lz)/Lj) * AGHQ(Q,Wphi,f_y_z2,z,r_prime,k,Alpha,gamma,data_j,ind.x0,ind.x1,lambda)
    

    S_L_a <- det(Lz) * AGHQ(Q,Wphi,lnf_a_fy_1,z,r_prime,k,Alpha,gamma,data_j,ind.x0,ind.x1,lambda) 
    H_L_a <- det(Lz) * AGHQ(Q,Wphi,lnf_a_fy_2,z,r_prime,k,Alpha,gamma,data_j,ind.x0,ind.x1,lambda) 
    S_L_rc <- det(Lz) * AGHQ(Q,Wphi,lnf_rc_fy_1,z,r_prime,k,Alpha,gamma,data_j,ind.x0,ind.x1,lambda)
    H_L_rc <- det(Lz) * AGHQ(Q,Wphi,lnf_rc_fy_2,z,r_prime,k,Alpha,gamma,data_j,ind.x0,ind.x1,lambda)
    H_L_a_rc <- det(Lz) * AGHQ(Q,Wphi,lnf_rc_a_fy_2,z,r_prime,k,Alpha,gamma,data_j,ind.x0,ind.x1,lambda)
    
    S_L_t <- det(Lz) * AGHQ(Q,Wphi,lnf_t_fy_1,z,r_prime,k,Alpha,gamma,data_j,ind.x0,ind.x1,lambda,tau=tau) #r(r+1)/2 x 1
    H_L_t <- det(Lz) * AGHQ(Q,Wphi,lnf_t_fy_2,z,r_prime,k,Alpha,gamma,data_j,ind.x0,ind.x1,lambda,tau=tau)
    H_L_at <- det(Lz) * AGHQ(Q,Wphi,lnf_t_a_fy_2,z,r_prime,k,Alpha,gamma,data_j,ind.x0,ind.x1,lambda,tau=tau)
    S_L_d <- det(Lz) * AGHQ(Q,Wphi,lnf_delta_fy_1,z,r_prime,k,Alpha,gamma,data_j,ind.x0,ind.x1,lambda,delta = delta) 
    H_L_d <- det(Lz) * AGHQ(Q,Wphi,lnf_delta_fy_2,z,r_prime,k,Alpha,gamma,data_j,ind.x0,ind.x1,lambda,delta = delta)
    H_L_ad <- det(Lz) * AGHQ(Q,Wphi,lnf_delta_a_fy_2,z,r_prime,k,Alpha,gamma,data_j,ind.x0,ind.x1,lambda, delta=delta)
    H_L_rc_t <- det(Lz) * AGHQ(Q,Wphi,lnf_t_rc_fy_2,z,r_prime,k,Alpha,gamma,data_j,ind.x0,ind.x1,lambda,tau = tau)
    H_L_rc_d <- det(Lz) * AGHQ(Q,Wphi,lnf_delta_rc_fy_2,z,r_prime,k,Alpha,gamma,data_j,ind.x0,ind.x1,lambda,delta = delta)
    H_L_t_d <- det(Lz) * AGHQ(Q,Wphi,lnf_delta_t_fy_2,z,r_prime,k,Alpha,gamma,data_j,ind.x0,ind.x1,lambda,tau = tau,delta = delta)
    
    S_L_lambda1 <- det(Lz) * AGHQ(Q,Wphi,lnf_lambda_fy_1,z,r_prime,k,Alpha,gamma,data_j,ind.x0,ind.x1,lambda,M1)
    H_L_lambda1 <- det(Lz) * AGHQ(Q,Wphi,lnf_lambda_fy_2,z,r_prime,k,Alpha,gamma,data_j,ind.x0,ind.x1,lambda,M1)
    H_L_a_lambda1 <- det(Lz) * AGHQ(Q,Wphi,lnf_a_lambda_fy_2,z,r_prime,k,Alpha,gamma,data_j,ind.x0,ind.x1,lambda,M1)
    H_L_rc_lambda1 <- det(Lz) * AGHQ(Q,Wphi,lnf_rc_lambda_fy_2,z,r_prime,k,Alpha,gamma,data_j,ind.x0,ind.x1,lambda,M1)
    H_L_t_lambda1 <- det(Lz) * AGHQ(Q,Wphi,lnf_t_lambda_fy_2,z,r_prime,k,Alpha,gamma,data_j,ind.x0,ind.x1,lambda, tau = tau, M1)
    H_L_d_lambda1 <- det(Lz) * AGHQ(Q,Wphi,lnf_delta_lambda_fy_2,z,r_prime,k,Alpha,gamma,data_j,ind.x0,ind.x1,lambda, delta, M1)
    
    if(nlambda==2){
      S_L_lambda2 <- det(Lz) * AGHQ(Q,Wphi,lnf_lambda_fy_1,z,r_prime,k,Alpha,gamma,data_j,ind.x0,ind.x1,lambda,M2)
      H_L_lambda2 <- det(Lz) * AGHQ(Q,Wphi,lnf_lambda_fy_2,z,r_prime,k,Alpha,gamma,data_j,ind.x0,ind.x1,lambda,M2)
      H_L_a_lambda2 <- det(Lz) * AGHQ(Q,Wphi,lnf_a_lambda_fy_2,z,r_prime,k,Alpha,gamma,data_j,ind.x0,ind.x1,lambda,M2)
      H_L_rc_lambda2 <- det(Lz) * AGHQ(Q,Wphi,lnf_rc_lambda_fy_2,z,r_prime,k,Alpha,gamma,data_j,ind.x0,ind.x1,lambda,M2)
      H_L_lambda1_lambda2 <- det(Lz) * AGHQ(Q,Wphi,lnf_lambda1_lambda2_fy_2,z,r_prime,k,Alpha,gamma,data_j,ind.x0,ind.x1,lambda,M12)
      H_L_t_lambda2 <- det(Lz) * AGHQ(Q,Wphi,lnf_t_lambda_fy_2,z,r_prime,k,Alpha,gamma,data_j,ind.x0,ind.x1,lambda, tau = tau, M2)
      H_L_d_lambda2 <- det(Lz) * AGHQ(Q,Wphi,lnf_delta_lambda_fy_2,z,r_prime,k,Alpha,gamma,data_j,ind.x0,ind.x1,lambda, delta, M2)
      
    }
    
    logLj <- log(Lj)

    
    Ha_j <- (H_L_a/Lj - S_L_a %*% t(S_L_a)/(Lj^2)) 
    Hrc_j <-  (H_L_rc/Lj - S_L_rc %*% t(S_L_rc)/(Lj^2)) 
    Harc_j <- (H_L_a_rc/Lj - S_L_rc %*% t(S_L_a)/(Lj^2)) 
    
    Ht_j <- (H_L_t/Lj - S_L_t %*% t(S_L_t)/(Lj^2))
    Hat_j <- (H_L_at/Lj - S_L_t %*% S_L_a/(Lj^2)) 
    Hrct_j <- (H_L_rc_t/Lj - t(S_L_rc %*% t(S_L_t))/(Lj^2)) 
    
    Hd_j <- (H_L_d/Lj - S_L_d %*% t(S_L_d)/(Lj^2)) 
    Had_j <- (H_L_ad/Lj - S_L_d %*% S_L_a/(Lj^2)) 
    Htd_j <- (H_L_t_d/Lj - t(S_L_t %*% t(S_L_d))/(Lj^2)) 
    Hrcd_j <- (H_L_rc_d/Lj - t(S_L_rc %*% t(S_L_d))/(Lj^2))
    
    Hlambda1_j <- (H_L_lambda1/Lj - S_L_lambda1 %*% t(S_L_lambda1)/(Lj^2))
    Halambda1_j <- (H_L_a_lambda1/Lj - S_L_lambda1 %*% t(S_L_a)/(Lj^2))
    Hrclambda1_j <- (H_L_rc_lambda1/Lj - S_L_lambda1 %*% t(S_L_rc)/(Lj^2))
    Htlambda1_j <- (H_L_t_lambda1/Lj - S_L_t %*% S_L_lambda1/(Lj^2))
    Hdlambda1_j <- (H_L_d_lambda1/Lj - S_L_d %*% S_L_lambda1/(Lj^2))
    
    Hlambda2_j <-0
    Halambda2_j <-0
    Hrclambda2_j <-0
    Hlambda1_2_j <-0
    Htlambda2_j <-0
    Hdlambda2_j <-0
    
    if(nlambda==2){

      Hlambda2_j <- (H_L_lambda2/Lj - S_L_lambda2 %*% t(S_L_lambda2)/(Lj^2))
      Halambda2_j <- (H_L_a_lambda2/Lj - S_L_lambda2 %*% t(S_L_a)/(Lj^2))
      Hrclambda2_j <- (H_L_rc_lambda2/Lj - S_L_lambda2 %*% t(S_L_rc)/(Lj^2))
      Hlambda1_2_j <- (H_L_lambda1_lambda2/Lj - S_L_lambda2 %*% t(S_L_lambda1)/(Lj^2))
      Htlambda2_j <- (H_L_t_lambda2/Lj - S_L_t %*% S_L_lambda2/(Lj^2))
      Hdlambda2_j <- (H_L_d_lambda2/Lj - S_L_d %*% S_L_lambda2/(Lj^2))
      
    }
    

      list2 <- list(E_u_j=E_u_j,E_u2=E_u2, logLj=logLj, 
                    Ha_j=Ha_j, Hrc_j=Hrc_j, Harc_j=Harc_j,
                    Ht_j=Ht_j, Hat_j=Hat_j, Hrct_j=Hrct_j,
                    Hd_j=Hd_j, Had_j=Had_j, Hrcd_j=Hrcd_j, Htd_j=Htd_j,
                    Hlambda1_j=Hlambda1_j, Halambda1_j=Halambda1_j, Hrclambda1_j=Hrclambda1_j, 
                    Htlambda1_j=Htlambda1_j,Hdlambda1_j=Hdlambda1_j,
                    Hlambda2_j=Hlambda2_j, Halambda2_j=Halambda2_j, Hrclambda2_j=Hrclambda2_j, Hlambda1_2_j=Hlambda1_2_j, 
                    Htlambda2_j=Htlambda2_j,Hdlambda2_j=Hdlambda2_j)
    
  }#end of loop for clinic
  logL <- 0
  tau.del.hat <-0
 
  Ha <- 0
  Harc <- 0
  Hrc <- 0
  Hat <- 0
  Ht <- 0
  Had <- 0
  Hrct <- 0
  Hrcd <- 0
  Hd <- 0
  Htd <- 0
  Hlambda1 <- 0
  Hlambda2 <- 0
  Halambda1 <- 0
  Halambda2 <- 0
  Hrclambda1 <- Hrclambda2 <- 0
  Hlambda1_2 <- 0
  Htlambda1 <- Htlambda2 <- Hdlambda1 <- Hdlambda2 <- 0
  
  for(j in 1:J) {
    tau.del.hat <- tau.del.hat+list1[[j]]$E_u2
    logL <- logL+list1[[j]]$logLj
    
    Ha <- Ha+list1[[j]]$Ha_j
    Hrc <- Hrc+list1[[j]]$Hrc_j
    Harc <- Harc+list1[[j]]$Harc_j
    
    Ht <- Ht + list1[[j]]$Ht_j
    Hat <- Hat + list1[[j]]$Hat_j
    Had <- Had + list1[[j]]$Had_j
    Hrct <- Hrct + list1[[j]]$Hrct_j
    Hrcd <- Hrcd + list1[[j]]$Hrcd_j
    Hd <- Hd + list1[[j]]$Hd_j
    Htd <- Htd + list1[[j]]$Htd_j
    
    Hlambda1 <- Hlambda1+list1[[j]]$Hlambda1_j
    Hlambda2 <- Hlambda2+list1[[j]]$Hlambda2_j
    Halambda1 <- Halambda1+list1[[j]]$Halambda1_j
    Halambda2 <- Halambda2+list1[[j]]$Halambda2_j
    Hrclambda1 <- Hrclambda1+list1[[j]]$Hrclambda1_j
    Hrclambda2 <- Hrclambda2+list1[[j]]$Hrclambda2_j
    Hlambda1_2 <- Hlambda1_2+list1[[j]]$Hlambda1_2_j
    
    Htlambda1 <- Htlambda1+list1[[j]]$Htlambda1_j
    Htlambda2 <- Htlambda2+list1[[j]]$Htlambda2_j
    Hdlambda1 <- Hdlambda1+list1[[j]]$Hdlambda1_j
    Hdlambda2 <- Hdlambda2+list1[[j]]$Hdlambda2_j
    
  }
  
  E_u <- list1[[1]]$E_u_j

  for(j in 2:J){
    E_u <- rbind(E_u, list1[[j]]$E_u_j)
  }
  
  print(paste0("loglike = ", logL))
  if(nlambda==1){
    H1 <- cbind(Ha,t(Harc),t(Halambda1),t(Hat),t(Had))
    H2 <- cbind(Harc,Hrc,t(Hrclambda1),t(Hrct),t(Hrcd))
    H3 <- cbind(Halambda1,Hrclambda1,Hlambda1,t(Htlambda1),t(Hdlambda1))
    H4 <- cbind(Hat,Hrct,Htlambda1,Ht,t(Htd))
    H5 <- cbind(Had,Hrcd,Hdlambda1,Htd,Hd)
    H <- rbind(H1,H2,H3,H4,H5)
  }
  
  if(nlambda==2){
    H1 <- cbind(Ha,t(Harc),t(Halambda1),t(Halambda2),t(Hat),t(Had))
    H2 <- cbind(Harc,Hrc,t(Hrclambda1),t(Hrclambda2),t(Hrct),t(Hrcd))
    H3 <- cbind(Halambda1,Hrclambda1,Hlambda1,t(Hlambda1_2),t(Htlambda1),t(Hdlambda1))
    H4 <- cbind(Halambda2,Hrclambda2,Hlambda1_2,Hlambda2, t(Htlambda2), t(Hdlambda2))
    H5 <- cbind(Hat,Hrct,Htlambda1,Htlambda2, Ht,t(Htd))
    H6 <- cbind(Had,Hrcd,Hdlambda1,Hdlambda2,Htd,Hd)
    H <- rbind(H1,H2,H3,H4,H5,H6)
  }
  
  print(H)
  print(paste0("symmetric H? ",isSymmetric(H, tol=1e-14)))
  se<-sqrt(diag(solve(-H)))
  
  pval <- c()
  for(i in 1:length(para0)){
    if(para0[i]>=0 & !is.nan(se[i])){
      pval[i] <- pnorm((para0[i]-0)/se[i],lower.tail = F) #para0 here is the converged estimates
    }
    if(para0[i]<0 & !is.nan(se[i])){
      pval[i] <- pnorm((para0[i]-0)/se[i],lower.tail = T) 
    }
  }
  
  #ac1-ac0
  cace <- para0[3]-para0[2]
  covc0c1 <- -solve(H)[2,3]
  cace.se <- sqrt(se[3]^2+se[2]^2-2*covc0c1)
  
  if(para0[3]-para0[2]>=0){
    cace.pval <- pnorm((para0[3]-para0[2])/cace.se,lower.tail = F) 
  } else {
    cace.pval <- pnorm((para0[3]-para0[2])/cace.se,lower.tail = T)
  }
  
  end_time <- Sys.time()
  
  time <- end_time-start_time
  print(time)
  
  return(list(par=para0,se=se,pval=pval,logL=logL,L2.norm=L2.norm[iter],time=time,cace=cace,cace.se=cace.se,cace.pval=cace.pval,niter=iter,E_u=E_u,Vu.b=Vu.b,H=H))
} #end of function

#The function returns par=estimates of parameters, se=standard errors of estimates, pval=p-values of wald test for estimates, logL=loglikelihood
#L2.norm=L2 norm, time=computation time, cace=mean cace, cace.se=standard error of CACE, cace.pval=p-value of wald test for cace, niter=# iterations, 
#E_u= site-level random effects of reduced dimension, Vu.b=vectorized covariance matrix of u and b from AGHQ (for calculating SE of site level CACE), H=Hessian matrix
