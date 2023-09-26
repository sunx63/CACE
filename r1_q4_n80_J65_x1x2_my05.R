#2023-5-18
#simulate population model when r=1 and estimate with model r=1
#dataset has 80 observations per site and 65 sites in total
#outcome variable Y is partially observed with 5% missing rate
#this is observed data, not the one with potential outcomes

rm(list=ls())

library(dplyr)
library(matrixcalc) #vec(), vech(), is.sysmetric.matrix
library(fastGHQuad) #integration
library(mvtnorm)
library(doParallel)
library(foreach)
library(lme4)

source("AGHQ.R")
source("f_y.R")
source("fun_der_lnf_alpha.R")
source("fun_2nd_deriv_lnf_a.R")
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

#simulation

miu.etaT=0.2
r=3
Q=4
S=1
r_prime=1
k=1
n <- 80
J=65
niter=500
tol <- 1e-04

alpha.true <- c(0.5,0.7,1.2,-0.5,1) #alpha_x1=1, alpha_x2=2
gamma.true <- c(1,-0.5, 0.5) #gamma_x1=0.5, gamma_x2=1
lambda.true <- c(1,0.7,1.1)
delta.true <- 0.3 #raise from 0.05
tau.true <- 0.5 # raise from 0.2

true <- c(alpha.true,gamma.true,lambda.true[2:3],tau.true, delta.true)

#functuion to similate data
sim_x1x2_2level <- function(seed, r, miu.etaT, n, J, alpha.true, gamma.true, lambda.true, tau.true, delta.true, S, missrate){
  set.seed(seed)
  alpha <- alpha.true
  gamma <- gamma.true
  N <- n*J
  
  #simulate T
  etaT=rnorm(J,miu.etaT,0.2) #
  piT=1/(1+exp(-etaT)) # P(T=1), constant within practice j
  summary(piT)   
  #   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  # 0.3540  0.4652  0.4982  0.5011  0.5339  0.6818
  id=T=rep(0,N)      # id, T, etaY, Y
  
  for (j in 1:J) {
    id[((j-1)*n+1):(j*n)] = j     # school ID
    T[((j-1)*n+1):(j*n)] = rbinom(n,1,piT[j]) # T~Bernoulli(piT[j])
  }

  x1 <- round(rnorm(N,1,1),4)
  x2 <- rbinom(N,1,0.65)
  
  #simulate C
  tau=delta.true  # var(delta)=tau 
  L=sqrt(tau)                       
  etaC=rep(0,N)  # etaC=r+b  compliance model
  
  Xc <- cbind(rep(1,N), x1, x2)
  
  b1 <- L * rnorm(J,0,1)
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
  #mean=0.97
  
  piY = 1/(1+exp(-etaY))
  
  Y=rep(0,N)  

  for(i in 1:N){
    Y[i] = rbinom(1,1,piY[i])
  } #generate outcome Y
  summary(Y)
  #mean=0.69. cor(Y,x2)=0.23, cor(Y,x1)=-0.21
  
  D=rep(0,N)  # set all D=0
  
  D[T*C[,2]==1]=1 # set D=1 of compliers assigned to T=1
  summary(D) #mean=0.39 
  summary(id)

  L1=as.data.frame(cbind(id, Y, T, D, C, x1, x2))   # L1 complete data=[id Y T D Cn Cc x1]
  names(L1)[5:6]=c("Cn", "Cc")    # Cn=1(C=n),Cc=1(C=c); for binary cases, C = Cc
  colnames(L1)[c(1,3)] = c("clinic","trt")
  L1$id <- rep(1:n,J)
  L1o=L1        

  L1o[T==0&D==0,5:6]=NA 
  
  #sample missing y
  miss.n <- missrate*N
  if(miss.n>0){
    miss.sample <- sample(1:N,miss.n,replace = F)
    L1o$M <- 1:N
    L1o$Y[which(L1o$M %in% miss.sample)] <- NA
    L1o <- L1o[,-length(L1o)]
  } else {
    L1o <- L1o
  }
  
  return(list(L1=L1,L1o=L1o)) #L1 is fully observed. L1o has missing data
}

#simulation
  print(paste0("ISIM= ",isim))
  L1o <- sim_x1x2_2level(seed = (isim), r=r_prime, miu.etaT, n, J, alpha.true, gamma.true, lambda.true, tau.true, delta.true, S=1, missrate=0.5)$L1o
  init.list <- set_init(L1o,r=3,side = 1,C~x1+x2+(1|clinic), Y~x1+x2+(1|clinic),"id","pmm",r_prime,col.clinic=1,col.trt=3,col.D=4,col.Y=2)

  init <- c(init.list$alpha.intercept.init, init.list$y.itt.est[2:length(init.list$y.itt.est)],
            init.list$C.init,
            init.list$lambda.init,
            init.list$tau.init,
            init.list$delta.init)
  #initial values should be in this order: alpha, gamma, lambda, tau, delta
  
  print(init)
  res <- cace_1side_parallel_lambda(r=r,k=k,Q=4,J=J,L1o,x0=c("x1","x2"),x1=c("x1","x2"),init = init, col.clinic=1,col.trt=3,col.D=4,col.Y=2,
                                    niter=niter,tol=tol,Share=1,r_prime=1)

  saveRDS(res,file = "r1my05_n80_J65_x1x2_isim.rds")

