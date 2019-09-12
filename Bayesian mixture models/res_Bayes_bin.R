## Script to collate results done by batches in HPC
## and obtain performace measures bias, CI coverage rate, etc

#install the necessary packages before loading them
library(compiler)
library(foreign)
#change directory to where your datasets and results are stored
source("perf.bayes.r")
perf.cmp <- cmpfun(perf.bayes)

beta.tc<-c(2,4)
N=c(1000, 200)
n.sims<-2500
#this specifies the length of the batch simulation done in HPC
everyout <- 25  


for (trt in 1:2){
  for (nmb in 1:2){
    for (scn in 1:2) {
      r.mu.cace <- r.se.cace<- r.lb.cace<- r.ub.cace <- r.med.cace  <-numeric(n.sims)
      
      
      for (i in 1:n.sims){
        
        # uncomment the line below to collate results for completely observed outcome
        #change the name by adding the suffix _wc0 for settings where there is a nonzero effect of complying in the controls
        
        for (i in seq(1,n.sims,everyout)){
          load(paste("binary_out/simulation_",beta.tc[trt],"/N",N[nmb],"/scenario_",scn,"/results/BayesUnobsX.missy.",i,"to",(i+everyout-1),".RData",sep=""))
          #load(paste("binary_out/simulation_",beta.tc[trt],"/N",N[nmb],"/scenario_",scn,"/results/BayesUnobsX.",i,"to",(i+everyout-1),".RData",sep=""))
          r.mu.cace[i:(i+everyout-1)]<-res$r.mu.cace[1:25]
          r.se.cace[i:(i+everyout-1)]<-res$r.se.cace[1:25]
          r.lb.cace[i:(i+everyout-1)]<-res$r.lb.cace[1:25]
          r.ub.cace[i:(i+everyout-1)]<-res$r.ub.cace[1:25]
          r.med.cace[i:(i+everyout-1)]<-res$r.med.cace[1:25]
          
          
        }
        
        results1<-as.data.frame(cbind(r.mu.cace,r.se.cace,r.lb.cace,r.ub.cace,r.med.cace))
        save(results1, file= paste("binary_out/simulation_",beta.tc[trt],"/N",N[nmb],"/scenario_",scn,"/results.Bayes.UX.missy.RData", sep="")) 
        #for binary outcome, when Late LOR=4, some datasets result in very large LOR even under the true princ strata
        #this are removed
        load(results,file=paste("binary_out/simulation_",beta.tc[trt],"/N",N[nmb],"/scenario_",scn,"/results_prin_cmp.Rdata",sep=""))  
        results1$r.cace<-results$r.cace
        result<-subset(results1,results1$r.cace<15)
        true.delta<-c(1.277708104, 2.748900682)
        perf.meas.Bayes<-perf.cmp(result$r.med.cace,result$r.se.cace,result$r.lb.cace,result$r.ub.cace,true.delta=true.delta[trt],n.sims=dim(result)[1])
        names(perf.meas.Bayes)<-c("coverage","bias","SE bias","rmse", "median ci width") 
        write.csv(perf.meas.Bayes, paste("Results/results_binout_sim",beta.tc[trt],"_N",N[nmb],"_sc",scn,"/perf_Bayes_missy.csv",sep=""))
        #remember to change the name of the file when collating fully obs outcome results
        #write.csv(perf.meas.Bayes, paste("Results/results_binout_sim",beta.tc[trt],"_N",N[nmb],"_sc",scn,"/perf_Bayes_fully.csv",sep=""))
        
        
      }
    }
  }
}
