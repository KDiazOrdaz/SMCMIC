setwd("~/Dropbox (Personal)/CACE MI")
#For trt=2, i.e. true LATE LOR =4 some simulated scenarios result in very high estimated LOR
#so we remove those from analyses 
N=c(1000, 200)
beta.tc<-c(2,4)
n.sims<-2500
r.cace<- r.se.cace  <-numeric(n.sims)

for (trt in 1:2){
  for (nmb in 1:2){
    for (scn in 1:2){
      
for (i in 1:n.sims){
        data<-read.dta(paste("binary_out/simulation_",beta.tc[trt],"/N",N[nmb],"/scenario_",scn,"/datasets/data_",i,".dta",sep=""))
        attach(data)
        rndcmp<-rnd*prin_cmp
        mod<-glm(y~ rndcmp + prin_cmp, family=binomial)
        beta <- coef(mod)[2]
      r.cace[i]<-beta
      r.se.cace[i]<-sqrt(diag(vcov(mod))[2])
      
      detach(data)

      results<-as.data.frame(cbind(r.cace,r.se.cace))
      
      save(results,file=paste("~/binary_out/simulation_",beta.tc[trt],"/N",N[nmb],"/scenario_",scn,"/results_prin_cmp.Rdata",sep=""))  
      
      
      }
}}}



###for settings where there is an effect of complying in the controls
r.cace<- r.se.cace  <-numeric(n.sims)
r.cace2<- r.se.cace2  <-numeric(n.sims)

for (trt in 1:2){
  for (nmb in 1:2){
    for (scn in 1:2){
      
      for (i in 1:n.sims){
        data<-read.dta(paste("binary_out/simulation_",beta.tc[trt],"/N",N[nmb],"/scenario_",scn,"/datasets/data_",i,".dta",sep=""))
        data$rndcmp<-data$rnd*data$prin_cmp
        mod<-glm(y~ rndcmp + prin_cmp, data=data,family=binomial)
        beta <- coef(mod)[2]
        r.cace[i]<-beta
        r.se.cace[i]<-sqrt(diag(vcov(mod))[2])
        
        results<-as.data.frame(cbind(r.cace,r.se.cace))
      
        save(results,file=paste("~/binary_out/simulation_",beta.tc[trt],"/N",N[nmb],"/scenario_",scn,"/results_prin_cmp.Rdata",sep=""))  
        
      }
    }}}

