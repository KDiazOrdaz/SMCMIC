library(compiler)
source("perf.bayes.r")
perf.cmp <- cmpfun(perf.bayes)
library(MASS) 
beta.tc<-c(2,4)
N=c(1000, 200)
n.sims<-2500
everyout <- 25  


for (trt in 1:2){
  for (nmb in 1:2){
    for (scn in 1:2) {
      r.mu.cace <- r.se.cace<- r.lb.cace<- r.ub.cace <- r.med.cace  <-numeric(n.sims)
      
      
      
      
      for (i in seq(1,n.sims,everyout)){
        #uncomment for collating results with fully observed outcome 
        #load(paste("cont_out/simulation_",beta.tc[trt],"/N",N[nmb],"/scenario_",scn,"/results/BayesUnobsX.fully.",i,"to",(i+everyout-1),".RData",sep=""))
        
        load(paste("cont_out/simulation_",beta.tc[trt],"/N",N[nmb],"/scenario_",scn,"/results/BayesUnobsX.missy,",i,"to",(i+everyout-1),".RData",sep=""))
        r.mu.cace[i:(i+everyout-1)]<-res$r.mu.cace[1:25]
        r.se.cace[i:(i+everyout-1)]<-res$r.se.cace[1:25]
        r.lb.cace[i:(i+everyout-1)]<-res$r.lb.cace[1:25]
        r.ub.cace[i:(i+everyout-1)]<-res$r.ub.cace[1:25]
        r.med.cace[i:(i+everyout-1)]<-res$r.med.cace[1:25]
              }
      
      results1<-as.data.frame(cbind(r.mu.cace,r.se.cace,r.lb.cace,r.ub.cace,r.med.cace))
      save(results1, file= paste("cont_out/simulation_",beta.tc[trt],"/N",N[nmb],"/scenario_",scn,"/results.Bayes.missy.RData", sep="")) 
      true.delta<-c(2, 4)
      perf.meas.Bayes<-perf.cmp(result$r.med.cace,result$r.se.cace,result$r.lb.cace,result$r.ub.cace,true.delta=true.delta[trt],n.sims=dim(results1)[1])
      names(perf.meas.Bayes)<-c("coverage","bias","SE bias","rmse", "median ci width") 
      write.csv(perf.meas.Bayes, paste("Results/results_contout_sim",beta.tc[trt],"_N",N[nmb],"_sc",scn,"/perf_Bayes_missy.csv",sep=""))
      #change for fully observed outcome
      #write.csv(perf.meas.Bayes, paste("Results/results_contout_sim",beta.tc[trt],"_N",N[nmb],"_sc",scn,"/perf_Bayes_fully.csv",sep=""))
      
      
    }
  }
}
