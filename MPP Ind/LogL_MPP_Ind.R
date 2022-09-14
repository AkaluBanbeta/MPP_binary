
#Log-likelihood function of modified power prior with independent power parameters
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

