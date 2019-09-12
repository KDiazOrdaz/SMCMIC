library(MASS) 
library(compiler)
source("performance.R")
perf.cmp <- cmpfun(performance)

beta.tc<-c(2,4)
N=c(1000, 200)
n.sims<-2500
  for (trt in 1:2){
    for (nmb in 1:2){
      for (scn in 1:2) {
        r.mu.cace <- r.se.cace  <-  numeric(n.sims)
            # collate results 
            for (i in 1:n.sims){
              try(load(paste("binary_out/simulation_",beta.tc[trt],"/N",N[nmb],"/scenario_",scn,"/results/SMC.results.unobservedX1.",i,".RData",sep="")), TRUE)
              
              r.mu.cace[i]<- results$r.mu.cace[1]
              r.se.cace[i]<- results$r.se.cace[1]
              results$r.mu.cace[1]<-NA
              results$r.se.cace[1]<-NA
	}
            
            results1<-as.data.frame(cbind(r.mu.cace,r.se.cace))
            save(results1, file= paste("binary_out/simulation_",beta.tc[trt],"/N",N[nmb],"/scenario_",scn,"/results.SMCFCS_UX_fully.RData", sep="")) 
            load(paste("binary_out/simulation_",beta.tc[trt],"/N",N[nmb],"/scenario_",scn,"/results_prin_cmp.Rdata",sep=""))  
            results1$r.cace<-results$r.cace
            result<-subset(results1,results1$r.cace<15)
            
            true.delta<-c(1.277708104, 2.748900682)
            
            perf<-perf.cmp(result$r.mu.cace,result$r.se.cace,true.delta=true.delta[2],n.sims=dim(results1)[1])
            names(perf)<-c("coverage","bias","SE bias","rmse", "median ci width") 
            write.csv(perf, paste("Results/results_binout_sim",beta.tc[trt],"_N",N[nmb],"_sc",scn,"/perf_SMCFCS_fully.csv",sep=""))
            
            
      }
    }
  }
