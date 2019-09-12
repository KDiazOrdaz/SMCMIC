library(MASS)
library("AER")
library(foreign)
library(boot)
library(mice)
library("mitools")
N=c(1000, 200)
beta.tc<-c(2,4)
n.sims<-2500
M<-10   #number of imputations

source("performance.R")


for (trt in 1:2){
  for (nmb in 1:2){
    for (scn in 1:2){
      r.mu.cace <- r.se.cace<-r.mu.cace.mi <- r.se.cace.mi <-r.mi.df <-numeric(n.sims)
        for (i in 1:n.sims){
          data<-read.dta(paste("cont_out/simulation_",beta.tc[trt],"/N",N[nmb],"/scenario_",scn,"/datasets/data_",i,".dta",sep=""))
          #need to create the treatment received var, 0 always in the control arm
          data$trt_rcv<-ifelse(data$prin_cmp,data$rnd,0) 
          #if two-way cross over then abs(1-rnd)
          
          ##analyse fully obs. outcome 
          data<-data[order(data$rnd),]   # order data by arm
          iv <- ivreg(data$y ~ data$trt_rcv | data$rnd)
          betaiv<-iv$coeff[2]
          se.iv<-sqrt(vcov(iv)[2,2])
          r.mu.cace[i]<-betaiv
          r.se.cace[i]<-se.iv
          ##analyse missing outcome, using multiple imputation by mice
          data.miss<-as.data.frame(cbind(data$y_miss, data$x2,data$rnd,data$trt_rcv))
          names(data.miss)<-c("y", "x2", "rnd","trt_rcv" )
          
          
          #MI using MICE in R 
          data<-data.miss[order(data.miss$rnd),]       # order data by arm
          #impute together 
          imp <-mice(data, m = M, method =c("norm","","",""),print=FALSE)
          impobj.miss <- imputationList(imp$impDatasets)
          
          iv <- with(imp, ivreg(y ~ trt_rcv  | rnd ))
          
          betaiv<-summary(pool(iv))[2,1]
          se.iv<-summary(pool(iv))[2,2]
        
          
        r.mu.cace.mi[i]<-betaiv
        r.se.cace.mi[i]<-se.iv
         } 
      
      results<-as.data.frame(cbind(r.mu.cace,r.se.cace))
      #save(results,file=paste("cont_out/simulation_",beta.tc[trt],"/N",N[nmb],"/scenario_",scn,"/results_IV2sls_unobsX1X2.Rdata",sep=""))  
      
      results1<-as.data.frame(cbind(r.mu.cace.mi,r.se.cace.mi,r.mi.df))
      #save(results1,file=paste("cont_out/simulation_",beta.tc[trt],"/N",N[nmb],"/scenario_",scn,"/results_IV2sls_ymiss.Rdata",sep=""))  
      
      
      true.delta<-beta.tc[trt]
      
      perf.meas<-performance(r.mu.cace, r.se.cace, true.delta,n.sims=n.sims)
      names(perf.meas)<-c("coverage","bias","MC error bias","rmse", "median ci width") 
      write.csv(perf.meas, paste("Results/results_contout_sim",beta.tc[trt],"_N",N[nmb],"_sc",scn,"/perf_TS_fully.csv",sep=""))

      perf.meas.mi<-performance(r.mu.cace.mi, r.se.cace.mi, true.delta,n.sims=n.sims)
      names(perf.meas.mi)<-c("coverage","bias","MC error bias","rmse", "median ci width") 
      write.csv(perf.meas.mi, paste("Results/results_contout_sim",beta.tc[trt],"_N",N[nmb],"_sc",scn,"/perf_TS_missy.csv",sep=""))
    } 
  }
}
