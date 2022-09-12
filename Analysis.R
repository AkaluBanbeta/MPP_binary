
# Analysing data based on the eight methods discussed in Banbeta et al (2019)

  # Data generation
  source("../Data_preparation.R")

  #Load functions for ccomputing the MPP
  source("../Functions_MPPs.R") 
  
  #1. MCMC sampling for the modified power prior with independent powers
  source("../Functions_MPPs.R") # function
  init_Ind=rep(0.3,H+2)
  result_MPP_Ind<-MCMCmetrop1R(LikelihoodMPP_Ind,theta.init=init_Ind,burnin =1*niter,mcmc =1*niter,verbose=0, thin = 1 ,V=0.0008*diag(length(init_Ind)),
                               success=success,H=H,n=n)
  result_MPP_Ind<-MCMCmetrop1R(LikelihoodMPP_Ind,theta.init=result_MPP_Ind[nrow(result_MPP_Ind),],  burnin = 1*niter, mcmc =0.8*niter,verbose=0, thin = 1 ,V=0.2*cov(result_MPP_Ind[trunc(0.8*niter):trunc(1*niter),]),
                               success=success,H=H,n=n)
  result_MPP_Ind<-MCMCmetrop1R(LikelihoodMPP_Ind,theta.init=result_MPP_Ind[nrow(result_MPP_Ind),],  burnin = 1*niter, mcmc =0.8*niter,verbose=0, thin = 1 ,V=0.4*cov(result_MPP_Ind[trunc(0.6*niter):trunc(0.8*niter),]),
                               success=success,H=H,n=n)
  result_MPP_Ind<-MCMCmetrop1R(LikelihoodMPP_Ind,theta.init=result_MPP_Ind[nrow(result_MPP_Ind),],  burnin = 1*niter, mcmc =1*niter,verbose=0, thin = 5 ,V=0.3*cov(result_MPP_Ind[trunc(0.6*niter):trunc(0.8*niter),]),
                               success=success,H=H,n=n)
  result_MPP_Ind<-as.mcmc(as.data.frame((result_MPP_Ind))) 
  varnames(result_MPP_Ind)[1]<-"theta"
  for(i in 1:H){
    varnames(result_MPP_Ind)[i+1]<-paste0("weight",i)
  }
  varnames(result_MPP_Ind)[H+2]<-"delta"
  summary(result_MPP_Ind)
  
  
  
  #2. MCMC sampling for for the modified power prior with depeendent power(DMPP)
  init_Dep=rep(0.3,H+4)
  result_MPP_Dep<-MCMCmetrop1R(LikelihoodMPP_Dep,theta.init=init_Dep,burnin =1*niter,mcmc =1*niter,verbose=0, thin = 1 ,V=0.0003*diag(length(init_Dep)),
                               success=success,H=H,n=n)
  result_MPP_Dep<-MCMCmetrop1R(LikelihoodMPP_Dep,theta.init=result_MPP_Dep[nrow(result_MPP_Dep),],  burnin = 1*niter, mcmc =0.8*niter,verbose=0, thin = 1 ,V=0.3*cov(result_MPP_Dep[trunc(0.8*niter):trunc(1*niter),]),
                               success=success,H=H,n=n)
  result_MPP_Dep<-MCMCmetrop1R(LikelihoodMPP_Dep,theta.init=result_MPP_Dep[nrow(result_MPP_Dep),],  burnin = 1*niter, mcmc =0.8*niter,verbose=0, thin = 1 ,V=0.22*cov(result_MPP_Dep[trunc(0.6*niter):trunc(0.8*niter),]),
                               success=success,H=H,n=n)
  result_MPP_Dep<-MCMCmetrop1R(LikelihoodMPP_Dep,theta.init=result_MPP_Dep[nrow(result_MPP_Dep),],  burnin = 1*niter, mcmc =1*niter,verbose=0, thin = 5 ,V=0.2*cov(result_MPP_Dep[trunc(0.6*niter):trunc(0.8*niter),]),
                               success=success,H=H,n=n)
  result_MPP_Dep<-as.mcmc(as.data.frame((result_MPP_Dep))) 
  varnames(result_MPP_Dep)[1]<-"p"
  for(i in 1:H){
    varnames(result_MPP_Dep)[i+1]<-paste0("weight",i)
  }
  varnames(result_MPP_Dep)[H+2]<-"a"
  varnames(result_MPP_Dep)[H+3]<-"b"
  varnames(result_MPP_Dep)[H+4]<-"delta"
  summary(result_MPP_Dep)
  
  
  #3. MCMC sampling for the Robust DMPP 1
  init_Robust=rep(0.3,H+4)
  result_MPP_Robust<-MCMCmetrop1R(LikelihoodMPP_Robust,theta.init=init_Robust,burnin =1*niter,mcmc =1*niter,verbose=0, thin = 1 ,V=0.0002*diag(length(init_Robust)),
                                  success=success,H=H,n=n)
  result_MPP_Robust<-MCMCmetrop1R(LikelihoodMPP_Robust,theta.init=result_MPP_Robust[nrow(result_MPP_Robust),],  burnin = 1*niter, mcmc =1*niter,verbose=0, thin = 1 ,V=0.3*cov(result_MPP_Robust[trunc(0.8*niter):trunc(1*niter),]),
                                  success=success,H=H,n=n)
  result_MPP_Robust<-MCMCmetrop1R(LikelihoodMPP_Robust,theta.init=result_MPP_Robust[nrow(result_MPP_Robust),],  burnin = 1*niter, mcmc =1*niter,verbose=0, thin = 1 ,V=0.25*cov(result_MPP_Robust[trunc(0.6*niter):trunc(0.8*niter),]),
                                  success=success,H=H,n=n)
  result_MPP_Robust<-MCMCmetrop1R(LikelihoodMPP_Robust,theta.init=result_MPP_Robust[nrow(result_MPP_Robust),],  burnin = 1*niter, mcmc =1*niter,verbose=0, thin = 5 ,V=0.15*cov(result_MPP_Robust[trunc(0.6*niter):trunc(0.8*niter),]),
                                  success=success,H=H,n=n)
  result_MPP_Robust<-as.mcmc(as.data.frame((result_MPP_Robust))) 
  varnames(result_MPP_Robust)[1]<-"p"
  for(i in 1:H){
    varnames(result_MPP_Robust)[i+1]<-paste0("weight",i)
  }
  varnames(result_MPP_Robust)[H+2]<-"a"
  varnames(result_MPP_Robust)[H+3]<-"b"
  varnames(result_MPP_Robust)[H+4]<-"delta"
  summary(result_MPP_Robust)
  
  
  
  #4. MCMC sampling for the Robust DMPP 2
  init_Robust_2=rep(0.3,H+4)
  result_MPP_Robust_2<-MCMCmetrop1R(LikelihoodMPP_Robust_2,theta.init=init_Robust_2,burnin =1*niter,mcmc =1*niter,verbose=0, thin = 1 ,V=0.0003*diag(length(init_Robust_2)),
                                    success=success,H=H,n=n)
  result_MPP_Robust_2<-MCMCmetrop1R(LikelihoodMPP_Robust_2,theta.init=result_MPP_Robust_2[nrow(result_MPP_Robust_2),],  burnin = 1*niter, mcmc =1*niter,verbose=0, thin = 1 ,V=0.35*cov(result_MPP_Robust_2[trunc(0.8*niter):trunc(1*niter),]),
                                    success=success,H=H,n=n)
  result_MPP_Robust_2<-MCMCmetrop1R(LikelihoodMPP_Robust_2,theta.init=result_MPP_Robust_2[nrow(result_MPP_Robust_2),],  burnin = 1*niter, mcmc =1*niter,verbose=0, thin = 1 ,V=0.35*cov(result_MPP_Robust_2[trunc(0.6*niter):trunc(0.8*niter),]),
                                    success=success,H=H,n=n)
  result_MPP_Robust_2<-MCMCmetrop1R(LikelihoodMPP_Robust_2,theta.init=result_MPP_Robust_2[nrow(result_MPP_Robust_2),],  burnin = 1*niter, mcmc =1*niter,verbose=0, thin = 5 ,V=0.2*cov(result_MPP_Robust_2[trunc(0.6*niter):trunc(0.8*niter),]),
                                    success=success,H=H,n=n)
  result_MPP_Robust_2<-as.mcmc(as.data.frame((result_MPP_Robust_2))) 
  varnames(result_MPP_Robust_2)[1]<-"p"
  for(i in 1:H){
    varnames(result_MPP_Robust_2)[i+1]<-paste0("weight",i)
  }
  varnames(result_MPP_Robust_2)[H+2]<-"a"
  varnames(result_MPP_Robust_2)[H+3]<-"b"
  varnames(result_MPP_Robust_2)[H+4]<-"delta"
  summary(result_MPP_Robust_2)
  
  
  #5. Bayesian analysis based on current data (no borrowing) 
  current_data=Data_list
  current_inits<- function() {list (theta_c=runif(1,0,0.5),delta=runif(1,0,0.49),.RNG.name="base::Mersenne-Twister",.RNG.seed=80000)}
  current_parameters.to.save<- c("delta","theta_c","theta_t")
  current_model<-jags.model(file="currentmodel.bug", data=current_data, inits=current_inits, n.chains=1, n.adapt=0.2*niter)
  update(current_model, n.iter=0.2*niter) # burn in
  result_current<-coda.samples(current_model, variable.names=current_parameters.to.save,thin = 5, n.iter=niter)
  varnames(result_current)[1]<-"delta"
  varnames(result_current)[2]<-"theta_c"
  varnames(result_current)[3]<-"theta_t"
  summary(result_current)
  
  
  
  #6. Bayesian pooled analysis
  data_pooled = Data_list
  inits_pooled <- function() {list (theta_c=runif(1,0,0.5),delta=runif(1,0,0.49),.RNG.name="base::Mersenne-Twister",.RNG.seed=80000)}
  pooled_parameters.to.save<- c("delta","theta_c","theta_t")
  pooled_model<-jags.model(file="pooledmodel.bug", data=data_pooled, inits=inits_pooled, n.chains=1, n.adapt=0.2*niter)
  update(pooled_model, n.iter=0.2*niter) # burn in
  result_pooled<-coda.samples(pooled_model, variable.names=pooled_parameters.to.save,thin = 5, n.iter=niter)
  varnames(result_pooled)[1]<-"delta"
  varnames(result_pooled)[2]<-"theta_c"
  varnames(result_pooled)[3]<-"theta_t"
  summary(result_pooled)
  
  
  #7. Bayesian meta-analytic-predictive (MAP) analysis
  MAP_data=Data_list
  inits_MAP =function() {list(delta = 0,tau = rexp(1, rate = 4),mu = rnorm(1, mean = -1, sd = 0.5),.RNG.name="base::Mersenne-Twister",.RNG.seed=80000)}
  MAP_parameters.to.save<- c("delta","tau","mu","theta_t")
  MAP_model<-jags.model(file="MAP_model.bug", data=MAP_data, inits=inits_MAP, n.chains=1, n.adapt=0.2*niter)
  update(MAP_model, n.iter=0.2*niter) # burn in
  result_MAP<-coda.samples(MAP_model, variable.names=MAP_parameters.to.save,thin = 5, n.iter=niter)
  varnames(result_MAP)[1]<-"delta"
  varnames(result_MAP)[2]<-"mu"
  varnames(result_MAP)[3]<-"tau"
  varnames(result_MAP)[4]<-"theta_t"
  summary(result_MAP)
  
  
  #8. Bayesian Robus MAP analysis
  Robust_MAP_data = Data_list
  inits_Robust_MAP =function() {list(delta = 0,tau = rexp(1, rate = 4),mu = rnorm(1, mean = -1, sd = 0.5),.RNG.name="base::Mersenne-Twister",.RNG.seed=80000)}
  Robust_MAP_parameters.to.save<- c("delta","tau","mu","theta_t")
  Robust_MAP_model<-jags.model(file="Robust_MAP_model.bug", data=Robust_MAP_data, inits=inits_Robust_MAP, n.chains=1, n.adapt=0.2*niter)
  update(Robust_MAP_model, n.iter=0.2*niter) # burn in
  result_Robust_MAP<-coda.samples(Robust_MAP_model, variable.names=Robust_MAP_parameters.to.save,thin = 5, n.iter=niter)
  varnames(result_Robust_MAP)[1]<-"delta"
  varnames(result_Robust_MAP)[2]<-"mu"
  varnames(result_Robust_MAP)[3]<-"tau"
  varnames(result_Robust_MAP)[4]<-"theta_t"
  summary(result_Robust_MAP)
  

  
  
