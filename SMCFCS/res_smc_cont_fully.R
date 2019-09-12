library(MASS) 
library(compiler)
source("performance.R")
perf.cmp <- cmpfun(performance)
#from the project directory
beta.tc<-c(2,4)
N=c(1000, 200)
n.sims<-2500
everyout <- 25  


  for (trt in 1:2){
    for (nmb in 1:2){
      for (scn in 1:2) {
        r.mu.cace <- r.se.cace <-numeric(n.sims)
        r.mu.cace.miss <- r.se.cace.miss <-numeric(n.sims)
               
        for (i in 1:n.sims){
      
            # collate results 
            for (i in seq(1,n.sims,everyout)){
              load(paste("cont_out/simulation_",beta.tc[trt],"/N",N[nmb],"/scenario_",scn,"/results/SMC.results.unobservedX1.",i,"to",(i+everyout-1),".RData",sep=""))

              r.mu.cace[i:(i+everyout-1)]<- results$r.mu.cace[1:25]
              r.se.cace[i:(i+everyout-1)]<- results$r.se.cace[1:25]
              
		

              
            }
            
            results1<-as.data.frame(cbind(r.mu.cace,r.se.cace))
            save(results1, file= paste("cont_out/simulation_",beta.tc[trt],"/N",N[nmb],"/scenario_",scn,"/results.SMCFCS.fully.RData", sep="")) 
            true.delta<-beta.tc
            
            
            perf<-perf.cmp(result$r.mu.cace,result$r.se.cace,true.delta=true.delta[trt],n.sims=dim(results1)[1])
            names(perf)<-c("coverage","bias","SE bias","rmse", "median ci width") 

            write.csv(perf, paste("Results/results_contout_sim",beta.tc[trt],"_N",N[nmb],"_sc",scn,"/perf_SMCFCS_fully.csv",sep=""))
            
         
      }
    }
  }
}
