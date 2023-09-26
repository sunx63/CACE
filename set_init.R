
#set_init() gives the initial values by fitting logit model by glmer()
#inputs:
#data: dataframe containing outcome, treatment assignment and receipt, cluster/site level variable 
#r: full dimension of random effects of Y model
#side: one sided or two sided noncompliance. default is one
#formula.c: formula for C model
#formula.y: formula for Y model
#subname: name of the cluster/site variable
#method: either "pmm" predictive mean matching or "prob" probability
#r_prime: reduced dimension of random effects of Y model
#col.clinic: column number of site-level variable
#col.trt: column number of treatment assignment
#col.D: column number of treatment receipt
#col.Y: column number of outcome

set_init <- function(data,r,side=1,formula.c,formula.y,subname,method,r_prime, col.clinic,col.trt,col.D,col.Y){
  cluster.ind <- which(colnames(data)==subname)
  colnames(data)[c(col.clinic, col.trt, col.D, col.Y)] <- c("clinic", "trt", "D", "Y") 
  if(side==1){
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
    }#end of loop
    P.c <- sum(data$C, na.rm = TRUE)/nrow(data[!is.na(data$C),])
    formC <- as.formula(formula.c)
    fit.c <- glmer(formC, data=data, family = binomial, nAGQ=20)
    C.init <- fit.c@beta
    delta.init <- fit.c@theta^2
    predicted <- predict(fit.c,type="response",newdata=data,allow.new.levels=T) #no random effect
    data$predicted <- predicted
    data$sub_id <- data[,cluster.ind]
    
    J <- dim(table(data$clinic))
    
    if(method=="pmm"){
      for(j in 1:J){
        subset <- data[data$clinic==j,]
        if(all(is.na(subset$C))) next
        for(i in 1:nrow(subset)){
          if(!is.na(subset$C[i])) next
          else {
            matchid <- subset[!is.na(subset$C),]$sub_id[which.min(abs(subset$predicted[i]-subset$predicted[!is.na(subset$C)]))]
            subset$C[i] <- subset$C[subset$sub_id==matchid]
            data$C[data$sub_id==subset$sub_id[i] & data$clinic==j] <- subset$C[i]
          }
        }#end of for loop i
      }#end of for loop j
    }#end of "pmm"
    
    if(method=="prob"){
      for(i in 1:nrow(data)){
        if(is.na(data$C[i])){
          data$C[i] <- data$predicted[i]
        }
      }
    }#end of "prob"
    if(r==1){
      fit.y <- glmer(Y~ -1 + I(1-C) + I((1-trt)*C) + I(trt*C)+(1|clinic), data=data, family=binomial)
    }
    if(r==2){
      fit.y <- glmer(Y~ -1 + I(1-C) + I((1-trt)*C) + I(trt*C)+(-1+I(1-C)+C|clinic), data=data, family=binomial)
    }
    if(r==3){
      fit.y <- glmer(Y~ -1 + I(1-C) + I((1-trt)*C) + I(trt*C)+(-1 + I(1-C) + I((1-trt)*C) + I(trt*C)|clinic), data=data, family=binomial)
    }

    alpha.init <- fit.y@beta
    
    formY <- as.formula(formula.y)
    fit.itt <- glmer(formY, data=data, family = binomial, nAGQ=20)

    if(r==1) {
      tau.init <- fit.y@theta^2
      lambda.init <- 0
    }
    if(r==2) {
      tau.vec <- as.data.frame(summary(fit.y)$varcor)[,4] #tau11, tau22, tau12
      tau.init <- matrix(c(tau.vec[1],tau.vec[3],tau.vec[3],tau.vec[2]),2,2)
      lambda.init <- 0
    }

    if(r==3) {
      if(r_prime==1){
        tau.vec <- as.data.frame(summary(fit.y)$varcor)[,4] #tau11, tau22, tau33, tau12, tau13, tau23
        tau.vec <- tau.vec[c(1,4,5,2,6,3)]
        tau.m <- VechToCovM(tau.vec, r)
        tau.init <- tau.m[1,1]
        lambda1.init <- tau.m[1,2]/tau.m[1,1]
        lambda2.init <- tau.m[1,3]/tau.m[1,1]
        lambda.init <- c(lambda1.init,lambda2.init)
      }
      if(r_prime==2){
        tau.vec <- as.data.frame(summary(fit.y)$varcor)[,4] #tau11, tau22, tau33, tau12, tau13, tau23
        tau.vec <- tau.vec[c(1,4,5,2,6,3)]
        tau.m <- VechToCovM(tau.vec, r)
        tau.init <- tau.m[1:2,1:2]
        lambda.init <- tau.m[1,3]/tau.m[1,1]
      }
      
      if(r_prime==3){
        tau.vec <- as.data.frame(summary(fit.y)$varcor)[,4] #tau11, tau22, tau33, tau12, tau13, tau23
        tau.vec <- tau.vec[c(1,4,5,2,6,3)]
        tau.init <- VechToCovM(tau.vec, r)
        lambda.init <- NA
      }
    } #end of if(r=3)

    y.itt.est <- fit.itt@beta
  }#end of if(side==1)

  return(list(alpha.intercept.init=alpha.init,tau.init=tau.init,C.init=C.init,delta.init=delta.init,lambda.init=lambda.init,
              estimatePC=P.c, y.itt.est=y.itt.est))
}#end of function

#set_init(L1o,r=r,side = 1,C~1+(1|clinic), Y~1+(1|clinic),"id","pmm",1)

#set_init(data,side = 1,C~age.c+first.colon.order+gender+(1|clinic), Y~trt*marry+trt*first.colon.order+(1|clinic))

#set_init(data, r=1, side = 1,C~x1+(1|clinic), Y~1+(1|clinic),"study_id")

