#Four functions for differet versions of the modified
#or normalized power prior (MPP) methods for binomial model


#1. MPP Ind: the MPP with Independent power parameters 
LikelihoodMPP_Ind <- function(parameters,success,H,n){
  theta=parameters[1]
  weight=parameters[2:(H+1)]    # weight is for the power (weight) parameter
  delta=parameters[H+2]         # delta for treatment effect
  theta2=delta+theta 
  Prior0=0
  if(any(weight>0 & weight<1) &theta>0 & theta<1 & theta2>0 & theta2<1){ 
    lik<-(sum(weight[1:H]*success[1:H],success[H+1]))*log(theta)+(sum(weight[1:H]*(n[1:H]-success[1:H]),(n[H+1]-success[H+1])))*log(1-theta)
    lik.1<-  success[H+2]*log(theta2)+(n[H+2]-success[H+2])*log(1-theta2) 
    Prior0<-Prior0+log(dbeta(theta,1,1 )) 
    Prior0<- Prior0+log(dnorm(delta,0,1000))
    Prior0<-Prior0+sum(log(dbeta(weight,1,1 )))
    ScalingConstantTemp <- lgamma(sum(weight[1:H]*success[1:H],1))+lgamma(sum(weight[1:H]*(n[1:H]-success[1:H]),1))-lgamma(sum(weight[1:H]*n[1:H],1,1))  
  }else{
    Prior0 <- -10^200
    lik.1<- -10^200
    lik <- -10^200
    ScalingConstantTemp<- 0
  }
  LikelihoodMPP_Ind<- lik+lik.1+Prior0-ScalingConstantTemp
}


# 2. The DMPP: the MPP with dependent power parameters 
LikelihoodMPP_Dep <- function(parameters,success,H,n){  
  theta=parameters[1]
  weight=parameters[2:(H+1)]
  a=parameters[H+2]
  b=parameters[H+3]
  delta=parameters[H+4]
  theta2=theta+delta
  Prior0=0
  mu=a/(a+b)
  sigmasq=(mu*(1-mu))/(a+b+1)
  if(any(weight >0 & weight<1) & theta>0 & theta<1 & theta2>0 & theta2<1 &
     mu>0 & mu <1 & sigmasq>0 & sigmasq<1 & a>0 & b>0){
    lik<-(sum(weight[1:H]*success[1:H])+success[H+1])*log(theta)+(sum(weight[1:H]*(n[1:H]-success[1:H]))+(n[H+1]-success[H+1]))*log(1-theta)
    lik.1<-  success[H+2]*log(theta2)+(n[H+2]-success[H+2])*log(1-theta2)
    Prior0<-Prior0+log(dbeta(theta,1,1))  
    Prior0<- Prior0+log(dnorm(delta,0,1000))
    Prior0<- Prior0+log(dunif(mu,0,1))
    Prior0<-Prior0+log(dgamma(1/sigmasq,0.01,0.01)/(sigmasq^2))
    Prior0<-Prior0+sum(log(dbeta(weight,a,b )))
    ScalingConstantTemp <- lgamma(sum(weight[1:H]*success[1:H])+1)+lgamma(sum(weight[1:H]*(n[1:H]-success[1:H]))+1)-lgamma(sum(weight[1:H]*n[1:H])+1+1)   
  }else{
    Prior0 <- -10^200
    lik <- -10^200
    lik.1<- -10^200
    ScalingConstantTemp<- 0
  }  
  LikelihoodMPP_Dep<- lik+lik.1+Prior0-ScalingConstantTemp
}


# 3. The Robust DMPP 1: robustification applied to each individual power parameter
LikelihoodMPP_Robust<-function(parameters,success,H,n){
  theta=parameters[1]
  weight=parameters[2:(H+1)]
  a=parameters[H+2]
  b=parameters[H+3]
  delta=parameters[H+4]
  theta2=theta+delta
  Prior0=0
  mu=a/(a+b)
  sigmasq=(mu*(1-mu))/(a+b+1)
  if(any(weight >0 & weight<1) & theta>0 & theta<1 & theta2>0 & theta2<1 &
     mu>0 & mu <1 & sigmasq>0 & sigmasq<1 & a>0 & b>0){
    lik<-(sum(weight[1:H]*success[1:H])+success[H+1])*log(theta)+(sum(weight[1:H]*(n[1:H]-success[1:H]))+(n[H+1]-success[H+1]))*log(1-theta)
    lik.1<-  success[H+2]*log(theta2)+(n[H+2]-success[H+2])*log(1-theta2)
    Prior0<-Prior0+log(dbeta(theta,1,1))  
    Prior0<- Prior0+log(dnorm(delta,0,1000))
    Prior0<- Prior0+log(dunif(mu,0,1))
    Prior0<-Prior0+log(dgamma(1/sigmasq,0.01,0.01)/(sigmasq^2))
    Prior0<-Prior0+sum(log(0.9*dbeta(weight,a,b )+0.1*dhalfnorm(weight,sd2theta(sqrt(sigmasq/6.25))))) 
    ScalingConstantTemp <- lgamma(sum(weight[1:H]*success[1:H])+1)+lgamma(sum(weight[1:H]*(n[1:H]-success[1:H]))+1)-lgamma(sum(weight[1:H]*n[1:H])+1+1) 
  }else{
    Prior0 <- -10^200
    lik <- -10^200
    lik.1<- -10^200
    ScalingConstantTemp<- 0
  }
  LikelihoodMPP_Robust<- lik+lik.1+Prior0-ScalingConstantTemp
}


#4. The Robust DMPP 2: robustification applied to the overall weights 
LikelihoodMPP_Robust_2<-function(parameters,success,H,n){
  theta=parameters[1]
  weight=parameters[2:(H+1)]
  a=parameters[H+2]
  b=parameters[H+3]
  delta=parameters[H+4]
  theta2=theta+delta
  Prior0=0
  mu=a/(a+b)
  sigmasq=(mu*(1-mu))/(a+b+1)
  if(any(weight >0 & weight<1) & theta>0 & theta<1 & theta2>0 & theta2<1 &
     mu>0 & mu <1 & sigmasq>0 & sigmasq<1 & a>0 & b>0){
    lik<-(sum(weight[1:H]*success[1:H])+success[H+1])*log(theta)+(sum(weight[1:H]*(n[1:H]-success[1:H]))+(n[H+1]-success[H+1]))*log(1-theta)
    lik.1<-  success[H+2]*log(theta2)+(n[H+2]-success[H+2])*log(1-theta2)
    Prior0<-Prior0+log(dbeta(theta,1,1))  
    Prior0<- Prior0+log(dnorm(delta,0,1000))
    Prior0<- Prior0+log(dunif(mu,0,1))
    Prior0<-Prior0+log(dgamma(1/sigmasq,0.01,0.01)/(sigmasq^2))
    Prior0<-Prior0+log (0.9*prod(dbeta(weight,a,b))+0.1*prod(dhalfnorm(weight,sd2theta(sqrt(sigmasq/6.25)))))
    ScalingConstantTemp <- lgamma(sum(weight[1:H]*success[1:H])+1)+lgamma(sum(weight[1:H]*(n[1:H]-success[1:H]))+1)-lgamma(sum(weight[1:H]*n[1:H])+1+1) 
  }else{
    Prior0 <- -10^200
    lik <- -10^200
    lik.1<- -10^200
    ScalingConstantTemp<- 0
  }  
  LikelihoodMPP_Robust_2<- lik+lik.1+Prior0-ScalingConstantTemp
}


