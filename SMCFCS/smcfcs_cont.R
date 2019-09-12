## script to perform Multiple imputation by SMC FCS of the compliance class
## on continuous outcomes
rm(list=ls())
#this runs in parallel each of the simulated datasets in each scenarios

# need to change the directory to where the dataset folder is
#"binary_out/simulation_",beta.tc[trt],"/N",N[nmb],"/scenario_",scn,
# for example 
## "CACE MI/cont_out/simulation_2/N1000/scenario_1/"
# get running counter from the system
taskID <- as.numeric(Sys.getenv("SGE_TASK_ID"))
print(taskID) 

# set specific values for the run 
everyout <- 25 # export results every 'everyout' runs 

require(foreign)
library(MASS)
library(mice)
library("mitools")
library(smcfcs)



n.sims <- 2500
count <- 0
numit=250
rjlimit=5000
m<-10
r.mu.cace <- r.se.cace <-numeric(n.sims)
r.mu.cace.miss <- r.se.cace.miss <-numeric(n.sims)


while (!file.exists(paste("results/SMC.results.unobservedX1.missy.",taskID,"to",(taskID+everyout-1),".RData",sep="")))
{
  for (i in taskID:(taskID+everyout-1)){
    count<-count+1
    cat("loop id= ",i,"\n")
  data<-read.dta(paste("datasets/data_",i,".dta",sep=""))

originaldata<-as.data.frame(cbind(data$y,data$rnd,data$cmp))
originaldata$cace<-data$rnd*data$cmp
names(originaldata)<-c("y", "rnd","cmp","cace" )


imp<-smcfcs(originaldata,smtype = "lm",
            smformula="y~  cmp + cace",
            method=c("","", "logreg","cmp*rnd"),
            predictorMatrix=NULL,
            m=m,
            numit=numit,
            rjlimit=rjlimit,
            noisy=FALSE) 

impobj <- imputationList(imp$impDatasets)
models <- with(impobj, lm(y~  cmp + cace))
mean.cace<-summary(MIcombine(models))[3,1]
se.cace<-summary(MIcombine(models))[3,2]

r.mu.cace[[count]]<-mean.cace
r.se.cace[[count]]<-se.cace

#missing
originaldata.miss<-as.data.frame(cbind(data$y_miss,data$x2,data$rnd,data$cmp))
originaldata.miss$cace<-data$rnd*data$cmp
names(originaldata.miss)<-c("y", "x2", "rnd","cmp","cace" )

imp.miss<-smcfcs(originaldata.miss,smtype = "lm",
                      smformula="y~  x2 + cmp + cace",
                      method=c("","","", "logreg","cmp*rnd"),
                      predictorMatrix=NULL,
                      m=m,
                      numit=numit,
                      rjlimit=rjlimit,
                      noisy=FALSE) 

  
impobj.miss <- imputationList(imp.miss$impDatasets)
models.miss <- with(impobj.miss, lm(y~  cmp + cace))

mean.cace.miss<-summary(MIcombine(models.miss))[3,1]
se.cace.miss<-summary(MIcombine(models.miss))[3,2]

r.mu.cace.miss[[count]]<-mean.cace.miss
r.se.cace.miss[[count]]<-se.cace.miss



} 
results<-as.data.frame(cbind(r.mu.cace,r.se.cace))
save(results,file=paste("results/SMC.results.unobservedX1.",taskID,"to",(taskID+everyout-1),".RData",sep=""))

results.miss<-as.data.frame(cbind(r.mu.cace.miss,r.se.cace.miss))
save(results.miss,file=paste("results/SMC.results.unobservedX1.missy.",taskID,"to",(taskID+everyout-1),".RData",sep=""))
}
