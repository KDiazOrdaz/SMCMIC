## script to perform Multiple imputation by SMC FCS of the compliance class
## on binary outcomes 
rm(list=ls())
#this runs in parallel each of the simulated datasets in each scenarios

# need to change the directory to where the dataset folder is
#"binary_out/simulation_",beta.tc[trt],"/N",N[nmb],"/scenario_",scn,
# for example 
## "CACE MI/binary_out/simulation_2/N1000/scenario_1/"


# get running counter from the system
taskID <- as.numeric(Sys.getenv("SGE_TASK_ID"))
print(taskID) 


# set specific values for the run 
everyout <- 0 # export results every 'everyout' runs 

require(foreign)
library(MASS)
library(mice)
library("mitools")
library(smcfcs)



n.sims <- 1
count <- 0
numit=250
rjlimit=5000
m<-10
r.mu.cace <- r.se.cace <-numeric(n.sims)
r.mu.cace.miss <- r.se.cace.miss <-numeric(n.sims)


while (!file.exists(paste("results/SMC.results.unobservedX1.missy.",taskID,"to",(taskID+everyout-1),".RData",sep="")))
  {
  #change below to load data in settings where there is an effect of complying in the controls.
  data<-read.dta(paste("datasets/data_",taskID,".dta",sep=""))
  
  originaldata<-as.data.frame(cbind(data$y,data$rnd,data$cmp))
  originaldata$cace<-data$rnd*data$cmp
  names(originaldata)<-c("y", "rnd","cmp","cace" )
  
  set.seed(10)
  imp<-smcfcs(originaldata,smtype = "logistic",
              smformula="y~  cmp + cace",
              method=c("","", "logreg","cmp*rnd"),
              predictorMatrix=NULL,
              m=m,
              numit=numit,
              rjlimit=rjlimit,
              noisy=FALSE) 
  
  impobj <- imputationList(imp$impDatasets)
  models <- with(impobj, glm(y~  cmp + cace, family=binomial()))
  mean.cace<-summary(MIcombine(models))[3,1]
  se.cace<-summary(MIcombine(models))[3,2]
  
  r.mu.cace[[taskID]]<-mean.cace
  r.se.cace[[taskID]]<-se.cace
  
  #missing
  originaldata.miss<-as.data.frame(cbind(data$y_miss,data$x2,data$rnd,data$cmp))
  originaldata.miss$cace<-data$rnd*data$cmp
  names(originaldata.miss)<-c("y", "x2", "rnd","cmp","cace" )
  
  imp.miss<-smcfcs(originaldata.miss,smtype = "logistic",
                   smformula="y~  x2 + cmp + cace",
                   method=c("","","", "logreg","cmp*rnd"),
                   predictorMatrix=NULL,
                   m=m,
                   numit=numit,
                   rjlimit=rjlimit,
                   noisy=FALSE) 
  
  
  impobj.miss <- imputationList(imp.miss$impDatasets)
  models.miss <- with(impobj.miss, glm(y~  cmp + cace, family=binomial()))
  
  mean.cace.miss<-summary(MIcombine(models.miss))[3,1]
  se.cace.miss<-summary(MIcombine(models.miss))[3,2]
  
  r.mu.cace.miss[[taskID]]<-mean.cace.miss
  r.se.cace.miss[[taskID]]<-se.cace.miss
  
  
  
  
  results<-as.data.frame(cbind(r.mu.cace,r.se.cace))
  save(results,file=paste("results/SMC.results.unobservedX1.",taskID,".RData",sep=""))
  
  results.miss<-as.data.frame(cbind(r.mu.cace.miss,r.se.cace.miss))
  save(results.miss,file=paste("results/SMC.results.unobservedX1.missy.",taskID,".RData",sep=""))
}
