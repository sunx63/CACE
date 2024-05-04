#2023-9-26
#simulate three-dimensioanl fully-observed Y, C partially missing.
#simulate J=170 sites and n=40 observations nested within 1 site.

rm(list=ls())
#setwd("/Users/elly/Documents/VCU/RA/NR/eassist/incomplete/script/missy")
library(dplyr)
library(matrixcalc) #vec(), vech(), is.sysmetric.matrix
library(fastGHQuad) #integration
library(mvtnorm)
library(doParallel)
library(foreach)
library(lme4)
library(reshape2)

source("AGHQ.R")
source("f_y.R")
source("deriv_lnf_alpha_1st.R")
source("deriv_lnf_alpha_2nd.R")
source("deriv_lnf_tau_1st.R")
source("deriv_lnf_delta_1st.R")
source("deriv_lnf_rc_1st.R")
source("fun_2nd_deriv_lnf_alpha_rc.R")
source("fun_2nd_deriv_lnf_rc.R")
source("deriv_lnf_lambda_1st.R")
source("deriv_lnf_lambda_2nd.R")
source("deriv_lnf_lambda_a_2nd.R")
source("deriv_lnf_lambda_rc_2nd.R")
source("deriv_lnf_lambda1_lambda2_2nd.R")
source("deriv_L_tau_2nd.R")
source("fun_2nd_L_delta.R")
source("VechToCovM.R")
source("MtoVech.R")
source("invMatrix.R")
source("set_init_noshare.R")
source("cace_1side_parallel.R")
registerDoParallel(cores=4)

#simulation
##fit intercept model 
miu.etaT=0.2
r=3
Q=4
r_prime=3 #three-dimensional Y model
k=1
n <- 40
J=170
niter=500
tol <- 1e-04

alpha.true <- c(0.5,0.7,1.2,-0.5,1) #alpha_x1=1, alpha_x2=2
gamma.true <- c(1,-0.5, 0.5) #gamma_x1=0.5, gamma_x2=1
lambda.true <- diag(c(1,1,1)) #r=3, no sharing, lambda is identity matrix
delta.true <- 0.3 
tau.true <- VechToCovM(c(0.3, 0.25,0.2,0.36, 0.4,0.6), r_prime) 

sim_x1x2_2level <- function(seed, r, miu.etaT, n, J, alpha.true, gamma.true, lambda.true, tau.true, delta.true){
  set.seed(seed)
  alpha <- alpha.true
  gamma <- gamma.true
  N <- n*J
  
  #simulate T
  etaT=rnorm(J,miu.etaT,0.2) #
  piT=1/(1+exp(-etaT)) # P(T=1), constant within site j
  summary(piT)   
  
  id=T=rep(0,N)      
  
  for (j in 1:J) {
    id[((j-1)*n+1):(j*n)] = j     # site ID
    T[((j-1)*n+1):(j*n)] = rbinom(n,1,piT[j]) # T~Bernoulli(piT[j])
  }
  
  x1 <- rnorm(N,1,1)
  x2 <- rbinom(N,1,0.65)
  
  #simulate C
  delta=delta.true 
 
  etaC=rep(0,N)  # etaC=r+b  compliance model
  
  Xc <- cbind(rep(1,N), x1, x2)
  
  b1 <-  rnorm(J,0,sqrt(delta))
  b2 <- list()
  for(j in 1:J){
    b2[[j]] <- rep(b1[j],n)
  }
  
  b <- unlist(b2)
  
  for(i in 1:N){
    etaC[i]=Xc[i,] %*% gamma +b[i]
  }
  summary(etaC)
  #mean=0.79
  
  pin=1/(1+exp(etaC))    # pin=P(C=n)
  pic=exp(etaC)/(1+exp(etaC)) 
  
  summary(pin)
  #mean=0.33
  summary(pic)
  #mean=0.67 
  
  C=matrix(0,N,2)        # C=[Cn Cc]=[1{C=n} 1{C=c}]
  for(i in 1:N){
    C[i,]=t(rmultinom(1,1,c(pin[i],pic[i])))
  } # C~binom(pic)
  summary(C)
  colnames(C) <- c("Cn","Cc")
  
  #simulate Y
  lambda <- lambda.true
  
  if(r==1){
    taubeta <- lambda %*% t(lambda) * tau.true
  } else {
    taubeta <- lambda %*% tau.true %*% t(lambda)
  }
  
  Lbeta=t(chol(taubeta))
  
  u=matrix(0,J,3)       #random effects of 3 compliance groups
  
  for (j in 1:J){
    u[j,] = t(Lbeta %*% rnorm(3,0,1))
  }
  
  u1 = u[rep(seq_len(nrow(u)), each = n),]
  
  Xy = cbind(x1,x2)
  
  beta=matrix(0,N,3) 
  
  for (i in 1:N){
    Aij <- cbind(diag(rep(1,3)), Xy[rep(i, 3),])
    beta[i,] <- Aij %*% alpha.true + u1[i,]
  }
  
  etaY=C[,1]*beta[,1]+C[,2]*((1-T)*beta[,2]+T*beta[,3])
  
  summary(etaY)
  
  piY = 1/(1+exp(-etaY))
  summary(piY)
  #mean=0.7076
  
  Y=rep(0,N)  
  
  for(i in 1:N){
    Y[i] = rbinom(1,1,piY[i])
  } #generate outcome Y
  
  summary(Y)
  #mean=0.69
  
  D=rep(0,N)  # set all D=0
  
  D[T*C[,2]==1]=1 # set D=1 of compliers assigned to T=1
  summary(D)
  #mean=0.38
  
  summary(id)
  
  #counterfactual
  T_cf <- 1-T
  Y_cf <- NA
  
  etaY_cf=C[,1]*beta[,1]+C[,2]*((1-T_cf)*beta[,2]+T_cf*beta[,3]) #need to add one more column of T for counter-factual scenario
  
  piY_cf = 1/(1+exp(-etaY_cf))
  summary(piY_cf)
  #mean=0.7049
  
  Y_cf=rep(0,N)  
  
  for(i in 1:N){
    Y_cf[i] = rbinom(1,1,piY_cf[i])
  } #generate outcome Y
  
  summary(Y)
  #mean=0.7041
  
  D_cf=rep(0,N)  # set all D=0
  
  D_cf[T_cf*C[,2]==1]=1 # set D=1 of compliers assigned to T=1
  summary(D_cf)
  #mean=0.29
  
  L1=as.data.frame(cbind(id, Y, T, D, C, Y_cf, T_cf, D_cf, x1, x2))   # L1 complete data=[id Y T D Cn Cc x1]
  
  colnames(L1)[c(1)] = c("clinic")
  
  L1$id <- rep(1:n,J)
  
  L1o=L1  #dataset with missing C  
  
  L1o[T==0&D==0,5:6]=NA 
  
  L2 <- L1[,c(1:3,6:8,10:12)]
  
  L2_long <- melt(L2,
                  id.vars=c("clinic", "Cc", "x1", "x2", "id"),
                  measure.vars=c("Y", "Y_cf"),
                  variable.name="y",
                  value.name="Y")
  L2_long <- L2_long[,-6]
  L2_long$trt <- c(T,T_cf)
  
  L3 <- L1[,c(1,3:4,6,8:12)]
  
  L3_long <- melt(L3,
                  id.vars=c("clinic", "Cc", "x1", "x2", "id"),
                  measure.vars=c("D", "D_cf"),
                  variable.name="d",
                  value.name="D")
  L2_long <- cbind(L2_long,L3_long[,7]) #dataset with potential outcomes, double the rows of L1
  colnames(L2_long)[ncol(L2_long)] <- "D"
  
  return(list(L1=L1,L1o=L1o, L2=L2_long))
} #end of function

  print(paste0("ISIM= ",isim))
  # isim is the loop id for simulations. isim goes from 1 to 500. Need a unix script to run multiple r simulations on computing clusters
  # if users just run a simulation, just manually change isim in the script to 1 or 2 or any numeric number. There are several locations of isim in the script

  L1o <- sim_x1x2_2level(seed = (isim), r=r_prime, miu.etaT, n, J, alpha.true, gamma.true, lambda.true, tau.true, delta.true)$L1o
  # "L1o" is the dataset with incomplete compliance

  init.list <- set_init_noshare(L1o,r=3,side = 1,C~x1+x2+(1|clinic), Y~x1+x2+(1|clinic),"id","pmm",col.clinic=1,col.trt=3,col.D=4,col.Y=2)
  # generate initial values
  # side is 1 for one-sided noncompliance
  # "id" is the subject level variable name. For this simulation, it is "id". Users need to adjust accordingly per their dataset.
  # pmm is predictive mean matching, please use as default.

  init <- c(init.list$alpha.intercept.init, init.list$y.itt.est[2:length(init.list$y.itt.est)],
            init.list$C.init,
            MtoVech(init.list$tau.init,r_prime),
            init.list$delta.init)
  #initial values should be in this order: alpha, gamma, tau, delta
  
  print(init)
  res <- cace_1side_parallel(r=r,k=k,Q=4,J=J,L1o,x0=c("x1","x2"),x1=c("x1","x2"),init = init, col.clinic=1,col.trt=3,col.D=4,col.Y=2,
                                    niter=niter,tol=tol,Share=0)

  saveRDS(res,file = "r3my0_n40_J170_x1x2_isim.rds")

