#2023-9-18
#potential outcomes

rm(list=ls())

library(lme4)
library(reshape2)
library(merDeriv) #vcov()

source("VechToCovM.R")
source("MtoVech.R")

#simulation

##fit random coefficient model 
miu.etaT=0.2
r=3
r_prime=1
k=1
n=40
J=170
niter=500
tol <- 1e-04

alpha.true <- c(0.5,0.7,1.2,-0.5,1) #alpha_x1=1, alpha_x2=2
gamma.true <- c(1,-0.5, 0.5) #gamma_x1=0.5, gamma_x2=1
lambda.true <- c(1,0.7,1.1)
delta.true <- 0.3 
tau.true <- 0.5 
converted.tau <- matrix(lambda.true) %*% tau.true %*% t(matrix(lambda.true))
converted.tau.vector <- MtoVech(converted.tau, 3)
converted.tau.vector
true <- c(alpha.true,converted.tau.vector,0.5,gamma.true,delta.true)


sim_x1x2_2level <- function(seed, r, miu.etaT, n, J, alpha.true, gamma.true, lambda.true, tau.true, delta.true){
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
   
  id=T=rep(0,N)      
  
  for (j in 1:J) {
    id[((j-1)*n+1):(j*n)] = j     # school ID
    T[((j-1)*n+1):(j*n)] = rbinom(n,1,piT[j]) # T~Bernoulli(piT[j])
  }
  
  x1 <- rnorm(N,1,1)
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



#simulation
nsim=500
#y model
npar.y <- 12
output_y <- matrix(NA,nrow = nsim,ncol=npar.y)
colnames(output_y) <- c("an","ac0","ac1","a1","a2","tau11","tau12","tau13","tau22","tau23","tau33","cace")
SE_y <- matrix(NA,nrow = nsim,ncol=npar.y)
colnames(SE_y) <- c("an","ac0","ac1","a1","a2","tau11","tau12","tau13","tau22","tau23","tau33","cace")
#c model
npar.c <- 4
output_c <- matrix(NA,nrow = nsim,ncol=npar.c)
colnames(output_c) <- c("r0","r1","r2","delta")
SE_c <- matrix(NA,nrow = nsim,ncol=npar.c)
colnames(SE_c) <- c("r0","r1","r2","delta")

options(max.print=100000)
options(nwarnings=10000)

        
time1 <- Sys.time()
for(isim in 1:nsim){
  print(paste0("sim= ", isim))
  L2 <- sim_x1x2_2level(seed = isim, r=r_prime, miu.etaT, n, J, alpha.true, gamma.true, lambda.true, tau.true, delta.true)$L2
  tryCatch({
  res <- glmer(Y~ -1 + I(1-Cc) + I((1-trt)*Cc) + I(trt*Cc)+x1+x2+(0+I(1-Cc) + I((1-trt)*Cc) + I(trt*Cc)|clinic), 
               data=L2, family=binomial)
  output_y[isim,c(1:5)] <- res@beta
  sd <- as.numeric(as.data.frame(VarCorr(res))[1:3,5]) #sd of random effects
  corr <- as.data.frame(VarCorr(res))[4:6,5]
  corrM <- VechToCovM(c(1,corr[1:2],1,corr[3],1),3) #correlation matrix of random effects
  tauM <- diag(sd) %*% corrM %*% diag(sd) #calculate tau estimates
  output_y[isim,6:11] <- MtoVech(tauM,3)
  output_y[isim,12] <- res@beta[3]-res@beta[2]
   
  SE_y[isim,1:5] <- summary(res)$coefficients[,2]
  SE_y[isim,6:11] <- sqrt(diag(vcov(res, full = TRUE)))[6:11] #SE for tau. the fixed effects se dont seem to match the previous ones
  SE_y[isim,12] <- sqrt(summary(res)$vcov[3,3]+summary(res)$vcov[2,2]-2*summary(res)$vcov[2,3])
  
  }, error=function(e){print(warnings())})
}
#
time2 <- Sys.time()
print(time2-time1)

for(isim in 1:nsim){
  print(paste0("sim= ", isim))
  L1 <- sim_x1x2_2level(seed = isim, r=r_prime, miu.etaT, n, J, alpha.true, gamma.true, lambda.true, tau.true, delta.true)$L1
  tryCatch({
    #c model
    res.c <- glmer(Cc~ x1+x2+(1|clinic), 
                   data=L1, family=binomial)
    output_c[isim,c(1:3)] <- res.c@beta
    output_c[isim,4] <- res.c@theta^2
    SE_c[isim,1:3] <- summary(res.c)$coefficients[,2]
    SE_c[isim,4] <- sqrt(diag(vcov(res.c, full = TRUE)))[4]
    
  }, error=function(e){print(warnings())})
}

#results
#setwd("/Users/elly/Documents/VCU/RA/NR/eassist/incomplete/sim/2023-9-18")
output <- cbind(output_y,output_c)
SE <- cbind(SE_y, SE_c)
write.csv(output,"output_potential_n40.csv", row.names = F)
write.csv(SE,"SE_potential_n40.csv", row.names = F)

est.potential <- round(apply(output,2,mean),3)
est.potential

ave.se <- round(apply(na.omit(SE),2,mean),3)
ave.se


cp.fun <- function(miuhat,se,miu){ #generate coverage probability
  upper <- miuhat+1.96*se
  lower <- miuhat-1.96*se
  coverage <- c()
  for(i in 1:length(upper)){
    if(!is.na(upper[i]) & !is.na(lower[i])){
      if(upper[i]>=miu & lower[i]<=miu){
        coverage[i] <- 1
      } else {
        coverage[i] <- 0
      }
    }
  }
  coverage.prob <- sum(na.omit(coverage))/sum(!is.na(coverage))
  return(coverage.prob)
}

#bias
bias <- est.potential-true 
varhat <- round(apply(output,2,var),3)

sdhat <- round(apply(output,2,sd),3)
sdhat

mse <- varhat + bias^2
relative.bias <- bias/true*100
relative.bias
cp <- c()

for(i in 1:ncol(output)){
  cp[i] <- cp.fun(output[,i],SE[,i],true[i])
}
r1.potential <- cbind(true,est.potential,ave.se, relative.bias, sdhat, mse, cp)
colnames(r1.potential) <- c("true","est", "se", "% bias", "sd", "mse", "CP")
write.csv(r1.potential, file = "r1_n40_potential_results.csv")




