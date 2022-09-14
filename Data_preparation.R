# This code generates binary data set with multiple trials

#Packages
library(parallel)
library(foreign)
library(Rlab)
library(stats4)
library(mcmc)
library(MCMCpack)
library(R2Cuba)
library(R2jags)
library(Hmisc)
library(gtools)
library(R.utils)
library(plyr)
library(dplyr)
library(R2WinBUGS) 
library("fdrtool")   

#Settings of simulation study
n_patients<- 100         #Number of patients in simulated data set
H<- 3                    #Number of historical trials
var_study_effect <- 0    #Variance of between-trials in log-odds historical studies.
pc <- 0.72
baseline_oddsratio <- log(0.72/0.28)    #Log baseline odds 
intervention_effect<- 0.13
niter=50000
#Treatment effect

sim=1
set.seed(sim+500)        #Use the same seed for each scenario

#Simulate the data set
trialnr<-rep(0,(H+2)*n_patients)
response<- rep(0,(H+2)*n_patients)
intervention<- rep(0,(H+2)*n_patients)


#Historical trials
for (i in 1:(H)){
  trial_effect<-rnorm(n=1,mean=0,sd=sqrt(var_study_effect)) + baseline_oddsratio
  trial_effect_indiv<- trial_effect
  for (j in 1:n_patients){    
    trialnr[(i-1)*n_patients+j]<-i
    response[(i-1)*n_patients+j]<-rbern(n=1,1/(1+exp(-trial_effect_indiv)))
  }  
}


#Current trial
i<- H+1
random_trial_effect_current_trial<-rnorm(n=1,mean=0,sd=sqrt(var_study_effect))
for (j in 1:n_patients){
  #Control group of current trial
  trial_effect<-random_trial_effect_current_trial+baseline_oddsratio
  trialnr[(i-1)*n_patients+j]<-H+1
  trial_effect_indiv<- trial_effect
  response[(i-1)*n_patients+j]<-rbern(n=1,1/(1+exp(-trial_effect_indiv)))
} 


#Intervention group of current trial
i<-H+2
for (j in 1:n_patients){
  trial_effect<-random_trial_effect_current_trial+log((intervention_effect+pc)/(1-(intervention_effect+pc)))
  intervention[(i-1)*n_patients+j]<-1
  trialnr[(i-1)*n_patients+j]<- H+1
  trial_effect_indiv<- trial_effect
  response[(i-1)*n_patients+j]<-rbern(n=1,1/(1+exp(-trial_effect_indiv)))
} 

dataset<-data.frame(trialnr=trialnr,response=factor(response,labels=c("alive","dead")),Intervention=intervention)
dataset<-data.frame(TrialArray=trialnr,response=response,InterventionArray=intervention)
dataset_long<-ddply(dataset,.variables=.(TrialArray,InterventionArray),summarise,success=sum(response))
dataset_long

success=dataset_long$success            # number of successes
n=rep(n_patients,H+2)       # number of patients per trial 

# This is the data we use in the analysis
Data_list=list(success=success,H=H,n=n)  

