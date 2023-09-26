
#run parellel jobs

#cace_1side_parallel() is the main function for estimating CACE with parallel computing with lambda=1. It is used when r=3 (full dimension)
#inputs
#r: full dimension of random effects of Y model. r=3 for one-sided noncompliance
#k: full dimension of random effects of C model. k=1 for one-sided noncompliance
#Q: number of abscissas for AGHQ
#J: number of clusters/sites
#x0: covariates controlled in Y model
#x1: covariates controlled in C model
#init: initial values for parameters
#data: dataframe containing outcome, treatment assignment and receipt, cluster/site level variable 
#col.clinic: column number of cluster/site level variable in the dataframe
#col.trt: column number of treatment assignment in the dataframe
#col.D: column number of treatment receipt in the dataframe
#col.Y: column number of outcome variable in the dataframe
#niter: maximum number of iteration. suggest 500
#tol: convergence criteria. suggest 10^-4
#Share: sharing of random effect. 
#When r=1, Share=1 means u'=un, Share=2 means u'=uc0 and Share=3 means u'=uc1
#When r=2, Share=1 means never takers and treatment compliers sharing one random effect, Share=2 means never takers and control compliers sharing one random effect, Share=3 means control compliers and trt compliers sharing one random effect

cace_1side_parallel <- function(r,k=1,Q,J,data,x0,x1,init,col.clinic,col.trt,col.D, col.Y, niter=500,tol,Share){
  start_time <- Sys.time()
  
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
  
  
  colnames(data)[c(col.clinic,col.trt, col.D, col.Y)] <- c("clinic", "trt", "D", "Y") 
  
  if(r==1){ #all 4 groups share a random effect
    lambda <- rep(1,3)
  }
  if(r==3){ 
    lambda <- diag(rep(1,r))
  }
  
  if(r==2){ 
    if(Share==0) lambda <- matrix(c(1,rep(0,3),1,1),nrow = 3,ncol = r,byrow = F) #compliers share a random effect
    if(Share==1) lambda <- matrix(c(rep(1,2),rep(0,3),1),nrow = 3,ncol = r,byrow = F) #never takers and control complier share a random effect
    if(Share==2) lambda <- matrix(c(1,0,1,0,1,0),nrow = 3,ncol = r,byrow = F) #never takers and trt compliers share a random effect
  }
  
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
  
  l.tau <- r*(r+1)/2
  l.d <- k*(k+1)/2
  
  theta <- matrix(NA,niter,(l.tau+l.alpha+l.gamma+l.d)) 
  #  colnames(theta) <- c('alpha_n','alpha_a','alpha_c0','alpha_c1','tou11','tou12','tou13','tou14','tou22','tou23','tou24','tou33','tou34','tou44',
  #gamma_a,gamma_c,delta11,delta12,delta22)
  theta[1,] <- init 
  
  if(r==1){
    Vu <- theta[1,(l.alpha+l.gamma+1):(l.tau+l.alpha+l.gamma)] 
  } else {
    Vu <- VechToCovM(theta[1,(l.alpha+l.gamma+1):(l.tau+l.alpha+l.gamma)], r)
  }
  
  if(k==1){
    Vb <- theta[1,((length(init)-l.d+1):length(init))]
  } else {
    Vb <- VechToCovM(theta[1,((length(init)-l.d+1):length(init))], k)
  }
  
  Vub <- bdiag(Vu, Vb) %>% as.matrix
  Vu.b <- do.call("rbind",replicate(J, MtoVech(Vub,(r+k)), simplify = FALSE))
  
  E_u <- do.call("rbind",replicate(J, rep(0,(r+k)), simplify = FALSE))
  
  E_u2 <- NA
  
  X <- fastGHQuad::gaussHermiteData(Q)
  
  #aghq matrix for u
  A <- matrix(NA, Q^(r+k),(r+k))
  
  for(i in 1:(r+k)){
    A[,i] <- rep(X$x,each=Q^((r+k)-i))
  }
  
  W <- matrix(NA, Q^(r+k),(r+k))
  for(i in 1:(r+k)){
    W[,i] <- rep(X$w,each=Q^((r+k)-i))
  }
  
  Ws <- W*exp(A^2)
  Ws <- apply(Ws,1,prod)
  
  L2.norm <- c() #l2 norm between theta^k and theta^(k+1)
  ll <- c() 

  for(iter in 1:niter){
    Alpha <- theta[iter,1:l.alpha]
    gamma <- theta[iter,(l.alpha+1):(l.alpha+l.gamma)]
    
    if(r==1){
      tau <- theta[iter,(l.alpha+l.gamma+1):(l.tau+l.alpha+l.gamma)]
    } else {
      tau <- VechToCovM(theta[iter,(l.alpha+l.gamma+1):(l.tau+l.alpha+l.gamma)],r)
    }
    
    if(k==1){
      delta <- theta[iter,(length(init)-l.d+1):length(init)]
    } else {
      delta <- VechToCovM(theta[iter,(length(init)-l.d+1):length(init)],k)
    }
    
    tau.delta <- bdiag(tau, delta) %>% as.matrix
    
    print(paste0("iteration = ", iter))
    
    # parallel jobs start here
   list1 <- foreach(j = 1:J,.verbose=T,.errorhandling='stop',.export=ls(.GlobalEnv),.packages=c('mvtnorm')) %dopar% { 
      #filter jth clinic data
      data_j <- data[data$clinic==j,]
      
      #change of variable
      V.z <- VechToCovM(Vu.b[j,],(k+r))
      
      Lz <- t(chol(2*V.z))
      z <- matrix(NA,Q^(r+k),(r+k))
      
      Wphi <- numeric(length=nrow(z))
      for(i in 1:nrow(z)){
        z[i,] <- Lz %*% A[i,]+E_u[j,]
        Wphi[i] <- Ws[i]*dmvnorm(z[i,], mean = rep(0, (r+k)), sigma = tau.delta, log = FALSE)
      }
      
      #AGHQ     
      Lj <- det(Lz)*AGHQ(Q,Wphi,f_y_incomplete,z,r,k,Alpha,gamma,data_j,ind.x0,ind.x1,lambda)

      E_u_j <- (det(Lz)/Lj) * AGHQ(Q,Wphi,f_y_z,z,r,k,Alpha,gamma,data_j,ind.x0,ind.x1,lambda)
      E_u2 <- (det(Lz)/Lj) * AGHQ(Q,Wphi,f_y_z2,z,r,k,Alpha,gamma,data_j,ind.x0,ind.x1,lambda)

      V <- E_u2 - E_u_j %*% t(E_u_j)
      
      Vub_j <- MtoVech(V,(r+k))
      
      S_L_a <- det(Lz) * AGHQ(Q,Wphi,lnf_a_fy_1,z,r,k,Alpha,gamma,data_j,ind.x0,ind.x1,lambda) 
      H_L_a <- det(Lz) * AGHQ(Q,Wphi,lnf_a_fy_2,z,r,k,Alpha,gamma,data_j,ind.x0,ind.x1,lambda) 
      S_L_rc <- det(Lz) * AGHQ(Q,Wphi,lnf_rc_fy_1,z,r,k,Alpha,gamma,data_j,ind.x0,ind.x1,lambda)
      H_L_rc <- det(Lz) * AGHQ(Q,Wphi,lnf_rc_fy_2,z,r,k,Alpha,gamma,data_j,ind.x0,ind.x1,lambda)
      H_L_a_rc <- det(Lz) * AGHQ(Q,Wphi,lnf_rc_a_fy_2,z,r,k,Alpha,gamma,data_j,ind.x0,ind.x1,lambda)
      
      logLj <- log(Lj)
      
      Sa_j <- S_L_a/Lj
      Src_j <- S_L_rc/Lj
      
      Ha_j <- (H_L_a/Lj - S_L_a %*% t(S_L_a)/(Lj^2)) #3*3
      Hrc_j <-  (H_L_rc/Lj - S_L_rc %*% t(S_L_rc)/(Lj^2)) #1
      Harc_j <- (H_L_a_rc/Lj - S_L_rc %*% t(S_L_a)/(Lj^2)) #1*3
      list2 <- list(E_u_j=E_u_j,E_u2=E_u2,Vub_j=Vub_j, logLj=logLj, Sa_j=Sa_j, Src_j=Src_j, Ha_j=Ha_j, Hrc_j=Hrc_j, Harc_j=Harc_j)
    }#end of loop for clinic
    
    tau.del.hat <- 0
    logL <- 0
    Sa <- 0
    Src <- 0
    Ha <- 0
    Hrc <- 0
    Harc <- 0
    
    for(j in 1:J) {
      tau.del.hat <- tau.del.hat+list1[[j]]$E_u2
      logL <- logL+list1[[j]]$logLj
      Sa <- Sa+list1[[j]]$Sa_j
      Src <- Src+list1[[j]]$Src_j
      Ha <- Ha+list1[[j]]$Ha_j
      Hrc <- Hrc+list1[[j]]$Hrc_j
      Harc <- Harc+list1[[j]]$Harc_j
    }
    
    tau.d.h <- tau.del.hat/J
    tau.h <- tau.d.h[1:r,1:r]
    d.h <- tau.d.h[(nrow(tau.d.h)-k+1):nrow(tau.d.h),(nrow(tau.d.h)-k+1):nrow(tau.d.h)]
    
    E_u <- list1[[1]]$E_u_j
    Vu.b <- list1[[1]]$Vub_j
 
    for(j in 2:J){
      E_u <- rbind(E_u, list1[[j]]$E_u_j)
      Vu.b <- rbind(Vu.b, list1[[j]]$Vub_j)
    }
    

    #log likelihood of J clinics
    print(paste0("loglike = ", logL))
    
    S <- c(Sa,Src)
    
    H1 <- cbind(Ha,t(Harc))
    H2 <- cbind(Harc,Hrc)
    H <- rbind(H1,H2)
    
    Hinv <- solve(H)
    
    theta[(iter+1),1:(l.alpha+l.gamma)] <- theta[iter,1:(l.alpha+l.gamma)] - as.vector(Hinv %*% S)
    
    if(r==1){
      theta[(iter+1),(length(Alpha)+length(gamma)+1):(l.tau+length(Alpha)+length(gamma))] <- tau.h
    } else {
      theta[(iter+1),(length(Alpha)+length(gamma)+1):(l.tau+length(Alpha)+length(gamma))] <- MtoVech(tau.h,r)
    }
    if(k==1){
      theta[(iter+1),((length(init)-l.d+1):length(init))] <- d.h
    } else {
      theta[(iter+1),((length(init)-l.d+1):length(init))] <- MtoVech(d.h,k)
    }
    
    print(paste0("theta= ",theta[(iter+1),]))
    
    if(iter>1){
      L2.norm[iter] <- sqrt(sum((theta[(iter),]-theta[(iter-1),])^2))
      print(paste0("L2.norm= ",L2.norm[iter]))
      if(L2.norm[iter] < tol) {break}
    }
    
    ll[iter] <- logL
    
  }

  para0 <- theta[(iter+1),]
  
  E_u2 <- NA
  
  Alpha <- para0[1:l.alpha]
  gamma <- para0[(length(Alpha)+1):(length(Alpha)+l.gamma)]
  
  if(r==1){
    tau <- para0[(length(c(Alpha,gamma))+1):(l.tau+length(Alpha)+length(gamma))]
  } else {
    tau <- VechToCovM(para0[(length(c(Alpha,gamma))+1):(l.tau+length(Alpha)+length(gamma))],r)
  }
  
  if(k==1){
    delta <- para0[(length(para0)-l.d+1):length(para0)]
  } else {
    delta <- VechToCovM(para0[(length(para0)-l.d+1):length(para0)],k)
  }
  
  tau.delta <- bdiag(tau, delta) %>% as.matrix
  

  list1 <- foreach(j = 1:J,.verbose=T,.errorhandling='stop',.export=ls(.GlobalEnv),.packages=c('mvtnorm')) %dopar% { 
    #filter jth clinic data
    data_j <- data[data$clinic==j,]
    
    #change of variable
    V.z <- VechToCovM(Vu.b[j,],(k+r))
    
    Lz <- t(chol(2*V.z))
    z <- matrix(NA,Q^(r+k),(r+k))
    
    Wphi <- numeric(length=nrow(z))
    for(i in 1:nrow(z)){
      z[i,] <- Lz %*% A[i,]+E_u[j,]
      Wphi[i] <- Ws[i]*dmvnorm(z[i,], mean = rep(0, (r+k)), sigma = tau.delta, log = FALSE)
    }
    
    #AGHQ     
    Lj <- det(Lz)*AGHQ(Q,Wphi,f_y_incomplete,z,r,k,Alpha,gamma,data_j,ind.x0,ind.x1,lambda)
    E_u_j <- (det(Lz)/Lj) * AGHQ(Q,Wphi,f_y_z,z,r,k,Alpha,gamma,data_j,ind.x0,ind.x1,lambda)
    E_u2 <- (det(Lz)/Lj) * AGHQ(Q,Wphi,f_y_z2,z,r,k,Alpha,gamma,data_j,ind.x0,ind.x1,lambda)
    V <- E_u2 - E_u_j %*% t(E_u_j)
    
    S_L_a <- det(Lz) * AGHQ(Q,Wphi,lnf_a_fy_1,z,r,k,Alpha,gamma,data_j,ind.x0,ind.x1,lambda) 
    H_L_a <- det(Lz) * AGHQ(Q,Wphi,lnf_a_fy_2,z,r,k,Alpha,gamma,data_j,ind.x0,ind.x1,lambda) 
    S_L_rc <- det(Lz) * AGHQ(Q,Wphi,lnf_rc_fy_1,z,r,k,Alpha,gamma,data_j,ind.x0,ind.x1,lambda)
    H_L_rc <- det(Lz) * AGHQ(Q,Wphi,lnf_rc_fy_2,z,r,k,Alpha,gamma,data_j,ind.x0,ind.x1,lambda)
    H_L_a_rc <- det(Lz) * AGHQ(Q,Wphi,lnf_rc_a_fy_2,z,r,k,Alpha,gamma,data_j,ind.x0,ind.x1,lambda)
    
    S_L_t <- det(Lz) * AGHQ(Q,Wphi,lnf_t_fy_1,z,r,k,Alpha,gamma,data_j,ind.x0,ind.x1,lambda,tau=tau) #r(r+1)/2 x 1
    H_L_t <- det(Lz) * AGHQ(Q,Wphi,lnf_t_fy_2,z,r,k,Alpha,gamma,data_j,ind.x0,ind.x1,lambda,tau=tau)
    H_L_at <- det(Lz) * AGHQ(Q,Wphi,lnf_t_a_fy_2,z,r,k,Alpha,gamma,data_j,ind.x0,ind.x1,lambda,tau=tau)
    S_L_d <- det(Lz) * AGHQ(Q,Wphi,lnf_delta_fy_1,z,r,k,Alpha,gamma,data_j,ind.x0,ind.x1,lambda,delta = delta) 
    H_L_d <- det(Lz) * AGHQ(Q,Wphi,lnf_delta_fy_2,z,r,k,Alpha,gamma,data_j,ind.x0,ind.x1,lambda,delta = delta)
    H_L_ad <- det(Lz) * AGHQ(Q,Wphi,lnf_delta_a_fy_2,z,r,k,Alpha,gamma,data_j,ind.x0,ind.x1,lambda, delta=delta)
    H_L_rc_t <- det(Lz) * AGHQ(Q,Wphi,lnf_t_rc_fy_2,z,r,k,Alpha,gamma,data_j,ind.x0,ind.x1,lambda,tau = tau)
    H_L_rc_d <- det(Lz) * AGHQ(Q,Wphi,lnf_delta_rc_fy_2,z,r,k,Alpha,gamma,data_j,ind.x0,ind.x1,lambda,delta = delta)
    H_L_t_d <- det(Lz) * AGHQ(Q,Wphi,lnf_delta_t_fy_2,z,r,k,Alpha,gamma,data_j,ind.x0,ind.x1,lambda,tau = tau,delta = delta)
    
    logLj <- log(Lj)
    
    Sa_j <- S_L_a/Lj
    Src_j <- S_L_rc/Lj
    St_j <- S_L_t/Lj
    Sd <- S_L_d/Lj 
    
    Ha_j <- (H_L_a/Lj - S_L_a %*% t(S_L_a)/(Lj^2)) #3*3
    Hrc_j <-  (H_L_rc/Lj - S_L_rc %*% t(S_L_rc)/(Lj^2)) #1
    Harc_j <- (H_L_a_rc/Lj - S_L_rc %*% t(S_L_a)/(Lj^2)) #1*3
    
    Ht_j <- (H_L_t/Lj - S_L_t %*% t(S_L_t)/(Lj^2))
    Hat_j <- (H_L_at/Lj - S_L_t %*% S_L_a/(Lj^2)) #1*3
    Hrct_j <- (H_L_rc_t/Lj - t(S_L_rc %*% t(S_L_t))/(Lj^2))
    
    Hd_j <- (H_L_d/Lj - S_L_d %*% t(S_L_d)/(Lj^2)) #1x1
    Had_j <- (H_L_ad/Lj - S_L_d %*% S_L_a/(Lj^2)) #1x3
    Hrcd_j <- (H_L_rc_d/Lj - t(S_L_rc %*% t(S_L_d))/(Lj^2)) #1x3
    Htd_j <- (H_L_t_d/Lj - t(S_L_t %*% t(S_L_d))/(Lj^2))
    
    list3 <- list(E_u_j=E_u_j, E_u2=E_u2, logLj=logLj, Ha_j=Ha_j, Hrc_j=Hrc_j, Harc_j=Harc_j, 
                  Ht_j=Ht_j, Hat_j=Hat_j, Hrct_j=Hrct_j, Hd_j=Hd_j, Had_j=Had_j, Hrcd_j=Hrcd_j, Htd_j=Htd_j)
  }#end of loop for clinic
  
  tau.del.hat <- 0
  logL <- 0

  Ha <- 0
  Hrc <- 0
  Harc <- 0
  Ht <- 0
  Hat <- 0
  Had <- 0
  Hrct <- 0
  Hrcd <- 0
  Hd <- 0
  Htd <- 0

  for(j in 1:J) {
    tau.del.hat <- tau.del.hat+list1[[j]]$E_u2
    logL <- logL + list1[[j]]$logLj
  
    Ha <- Ha + list1[[j]]$Ha_j
    Hrc <- Hrc + list1[[j]]$Hrc_j
    Harc <- Harc + list1[[j]]$Harc_j
    Ht <- Ht + list1[[j]]$Ht_j
    Hat <- Hat + list1[[j]]$Hat_j
    Had <- Had + list1[[j]]$Had_j
    Hrct <- Hrct + list1[[j]]$Hrct_j
    Hrcd <- Hrcd + list1[[j]]$Hrcd_j
    Hd <- Hd + list1[[j]]$Hd_j
    Htd <- Htd + list1[[j]]$Htd_j
  }
  
  tau.d.h <- tau.del.hat/J
  tau.h <- tau.d.h[1:r,1:r]
  d.h <- tau.d.h[(nrow(tau.d.h)-k+1):nrow(tau.d.h),(nrow(tau.d.h)-k+1):nrow(tau.d.h)]
  
  E_u <- list1[[1]]$E_u_j

  for(j in 2:J){
    E_u <- rbind(E_u, list1[[j]]$E_u_j)
  }
  
  print(paste0("loglike = ", logL))
  
  H1 <- cbind(Ha,t(Harc),t(Hat),t(Had))
  H2 <- cbind(Harc,Hrc,t(Hrct),t(Hrcd))
  H3 <- cbind(Hat,Hrct,Ht,t(Htd))
  H4 <- cbind(Had,Hrcd,Htd,Hd)
  H <- rbind(H1,H2,H3,H4)
  if(!is.symmetric.matrix(H)) print("H is not symmetric")
  print(H)
  
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
  covc0c1 <- solve(-H)[2,3]
  cace.se <- sqrt(se[3]^2+se[2]^2-2*covc0c1)
  
  if(para0[3]-para0[2]>=0){
    cace.pval <- pnorm((para0[3]-para0[2])/cace.se,lower.tail = F) #
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