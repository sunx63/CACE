#2022-11-29
#simulate e-assist data, J=170, number of patients nested within a physician varies from 1 to 44 according to the real study. 
#The patient distribution across physicians is stored in physician_order.rds. 

rm(list=ls())
#setwd("/Users/elly/Documents/VCU/RA/NR/eassist/incomplete/script/lambda_NR")
library(dplyr)
library(matrixcalc) #vec(), vech(), is.sysmetric.matrix
library(fastGHQuad) #integration
library(mvtnorm)
library(doParallel)
library(foreach)
library(lme4)

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
source("set_init.R")
source("cace_1side_parallel_lambda.R")

registerDoParallel(cores=40)
# users can define how many cores available according to their machines. 
# if the available cores < 40 or # specified in the above code, the machine will allocate what's avalable to the job.

#simulation
miu.etaT=0.2 # P(T=1) \approx 0.5
r=3 # for one-sided compliance, r=3
Q=4

r_prime=1 # reduced dimension of random effects in Y model
k=1 # for one-sided compliance, k=1
physician_order=readRDS("physician_order.rds") 
#e-assist study sample size file, required for e-assist study simulation
J=170 #sites
niter=500
tol <- 1e-04

alpha.true <- c(0.5,0.7,1.2,-0.5,1) #alpha_x1=1, alpha_x2=2
gamma.true <- c(1,-0.5, 0.5) #gamma_x1=0.5, gamma_x2=1
lambda.true <- c(1,0.7,1.1)
delta.true <- 0.3 
tau.true <- 0.5 

sim_x1x2_2level <- function(seed, r, miu.etaT, physician_order, J, alpha.true, gamma.true, lambda.true, tau.true, delta.true){
  set.seed(seed)
  alpha <- alpha.true
  gamma <- gamma.true
  N <- sum(physician_order)
  
  #simulate T
  etaT=rnorm(J,miu.etaT,0.2) #
  piT=1/(1+exp(-etaT)) # P(T=1), constant within practice j
  summary(piT)   
  #   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  # 0.3540  0.4652  0.4982  0.5011  0.5339  0.6818
  
  piT1 <- list()
  for(j in 1:J){
    piT1[[j]] <- rep(piT[j],physician_order[j])
  }
  piT1 <- unlist(piT1)
  
  T <- numeric(length=N)
  
  for(i in 1:N){
    T[i] <- rbinom(1,1,piT1[i])
  }

  id1 <- c(1:J)  
  id <- list()
  for(j in 1:J){
    id[[j]] <- rep(id1[j],physician_order[j])
  }
  id <- unlist(id)
  
  x1 <- rnorm(N,1,1)
  x2 <- rbinom(N,1,0.65)
  
  #simulate C
  delta=delta.true 
  
  etaC=rep(0,N)  # etaC=r+b  compliance model
  
  Xc <- cbind(rep(1,N), x1, x2)
  
  b1 <-  rnorm(J,0,sqrt(delta)) #delta is one-dimensional
  b2 <- list()
  for(j in 1:J){
    b2[[j]] <- rep(b1[j],physician_order[j])
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
  
  u1 = u[rep(seq_len(nrow(u)), physician_order),]

  Xy = cbind(x1,x2)
  
  beta=matrix(0,N,3) 
  
  for (i in 1:N){
    Aij <- cbind(diag(rep(1,3)), Xy[rep(i, 3),])
    beta[i,] <- Aij %*% alpha.true + u1[i,]
  }
  
  etaY=C[,1]*beta[,1]+C[,2]*((1-T)*beta[,2]+T*beta[,3])
  
  summary(etaY)
  #mean=0.97
  
  piY = 1/(1+exp(-etaY))
  
  Y=rep(0,N)  

  for(i in 1:N){
    Y[i] = rbinom(1,1,piY[i])
  } #generate outcome Y
  summary(Y)
  #mean=0.69. cor(Y,x2)=0.18, cor(Y,x1)=-0.19
  
  D=rep(0,N)  # set all D=0
  
  D[T*C[,2]==1]=1 # set D=1 of compliers assigned to T=1
  summary(D) #mean=0.39 
  summary(id)

  L1=as.data.frame(cbind(id, Y, T, D, C, x1, x2))   # L1 complete data=[id Y T D Cn Cc x1]
  names(L1)[5:6]=c("Cn", "Cc")    
  colnames(L1)[c(1,3)] = c("clinic","trt")

  ID <- list()
  for(j in 1:J){
    ID[[j]] <- c(1:physician_order[j])
  }
  ID <- unlist(ID)
  L1 <- cbind(L1, ID)
  
  L1o=L1        
  # L1o observed data=[id Y T D Cn Cc]
  
  L1o[T==0&D==0,5:6]=NA 
  return(list(L1=L1,L1o=L1o))
}


  print(paste0("ISIM= ",isim)) 
  # isim is the loop id for simulations. isim goes from 1 to 500. Need a unix script to run multiple r simulations on computing clusters
  # if users just run a simulation, just manually change isim in the script to 1 or 2 or any numeric number. There are several locations of isim in the script

  L1o <- sim_x1x2_2level(seed = (isim), r=r_prime, miu.etaT, physician_order, J, alpha.true, gamma.true, lambda.true, tau.true, delta.true)$L1o
  #L1o is the dataset of a multisite trial with e-assist sample size with observed/incomplete compliance

  # Generate initial values for the following algorithm
  init.list <- set_init(L1o,r=3,side = 1,C~x1+x2+(1|clinic), Y~x1+x2+(1|clinic),"ID","pmm",r_prime,col.clinic=1,col.trt=3,col.D=4,col.Y=2) 
  # this is to generate the initial values
  # side is 1 for one-sided noncompliance
  # ID is the subject level variable name. For this simulation, it is "ID". Users need to adjust accordingly per their dataset.
  # pmm is predictive mean matching, please use as default.
  # col.clinic, col.trt, col.D, col.Y are the column numbers of cluster ID, treatment assignment, treatment receipt/participation, outcome in the dataset

  
  init <- c(init.list$alpha.intercept.init, #initial values for an, ac0, ac1 in Y model, treatment effects of three compliance groups
            init.list$y.itt.est[2:length(init.list$y.itt.est)], #initial values for a1 and a2, coefficients for covariates x1 and x2 in Y model
            init.list$C.init, #initial values for r0, r1, r2, coefficients for intercept, x1 and x2 in C model
            init.list$lambda.init, #initial values for factor loading matrix
            init.list$tau.init, #initial values for cluster level variance of Y model
            init.list$delta.init) #initial values for cluster level variance of C model
  
  #Initial values should be in this order: alpha, gamma, lambda, tau, delta. Changing order will crash the program
  
  print(init)

  # estimation
  res <- cace_1side_parallel_lambda(r=r,k=k,Q=4,J=J,data=L1o,x0=c("x1","x2"),x1=c("x1","x2"),init = init, col.clinic=1,col.trt=3,col.D=4,col.Y=2, 
                                    niter=niter,tol=tol,Share=1,r_prime=1)
  # r: full dimension of random effects in Y model, numeric. r=3 for one-sided noncompliance.
  # k: full dimension of random effects in C model, numeric. k=1 for one-sided noncompliance.
  # Q: number of abscissas for AGHQ, numeric
  # J: number of sites, numeric
  # data: dataset
  # x0: covariates controlled in Y model, character
  # x1: covariates controlled in C model, character
  # init: intial values, numeric
  # col.clinic: column# of site/cluster level variable, numeric
  # col.trt: column# of treatment assignment, numeric
  # col.D: column# of treatment receipt, numeric
  # col.Y: column# of outcome, numeric
  # niter: # of iteration. program will stop if reaches niter specified, numeric. 300 is enough for e-assist study.
  # tol: tolerance for L2 norm, numeric. default is 10^-4
  # Share: sharing method, numeric. when r=1, Share=1, 2, 3 means never taker's, control complier's and treatment complier's random effect is shared among three groups respectively
  # r_prime: reduced dimension of random effects after sharing in Y model, numeric. when r=1, r_prime is 1, indicating one random effect is shared in Y model.

saveRDS(res,file = "r1_eassist_x1x2_isim.rds")

