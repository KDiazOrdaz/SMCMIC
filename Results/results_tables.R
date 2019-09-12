rm(list=ls())
out <- c("cont","bin")
ymiss <- c("fully","missy")
Method<-c("Bayes", "TS", "SMCFCS")
beta.tc<-c(2,4)
N=c(1000, 200)
NonCompliance<- c(30,50)
results<-list()
for (outc in 1:2) {
  results[[outc]]<-list()
  for (meth in 1:3){
  for (ymis in 1:2){
      for (trt in 1:2){
         for (nmb in 1:2){
           for (scn in 1:2) {
             
        df.t <- read.csv(paste("Results/results_",out[outc],"out_sim",beta.tc[trt],"_N",N[nmb],"_sc",scn,"/perf_",Method[meth],"_",ymiss[ymis],".csv",sep=""))
        df <- as.data.frame(t(df.t[,-1]))
        colnames(df)<-c("coverage",	"bias",	"mc.error",	"RMSE",	"CI width")
        
        df$Method<-Method[meth]
        df$MissingY<-ymiss[ymis]
        df$Effect.size<-beta.tc[trt]	
        df$Sample.Size<-N[nmb]	
        df$Noncompliance<-NonCompliance[scn]	
       
        results[[outc]]<-rbind(results[[outc]],df)
        
        }
    }
   }
  }
 }
}

#extra table for the binary outcomes with compliance class effect on the controls
for (outc in 2:2) {
  results[[3]]<-list()
  for (meth in 1:3){
    for (ymis in 1:2){
      for (trt in 1:2){
        for (nmb in 1:2){
          for (scn in 1:2) {
            
            df.t <- read.csv(paste("Results/results_",out[outc],"out_sim",beta.tc[trt],"_N",N[nmb],"_sc",scn,"/perf_",Method[meth],"_",ymiss[ymis],"_wc0.csv",sep=""))
            df <- as.data.frame(t(df.t[,-1]))
            colnames(df)<-c("coverage",	"bias",	"mc.error",	"RMSE",	"CI width")
            
            df$Method<-Method[meth]
            df$MissingY<-ymiss[ymis]
            df$Effect.size<-beta.tc[trt]	
            df$Sample.Size<-N[nmb]	
            df$NonCompliance<-NonCompliance[scn]	
            
            results[[3]]<-rbind(results[[3]],df)
          }
        }
      }
    }
  }
}
write.csv(as.data.frame(results[[1]]), file="Results/cont_results.csv",row.names=FALSE)
write.csv(as.data.frame(results[[2]]), file="Results/bin_results.csv",row.names=FALSE)
write.csv(as.data.frame(results[[3]]), file="Results/bin_results_wc0.csv",row.names=FALSE)
