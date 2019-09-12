#### two stage residual inclusion for fully observed binary outcome  

library(MASS)
library("AER")
library(foreign)
library(boot)
library(systemfit)
N=c(1000, 200)
beta.tc<-c(2,4)
n.sims<-2500
true.delta<-c(1.277708104, 2.748900682)

source("performance.R")


for (trt in 1:2){
  for (nmb in 1:2){
    for (scn in 1:2){
      r.mu.cace <- r.se.cace  <-numeric(n.sims)
      for (i in 1:n.sims){
        #change the nameof file by adding suffix _wc0 for the setting where there is a nonzero effect of complying in the controls
        data<-read.dta(paste("binary_out/simulation_",beta.tc[trt],"/N",N[nmb],"/scenario_",scn,"/datasets/data_",i,".dta",sep=""))
        #need to create the treatment received var, 0 always in the control arm
        data$trt_rcv<-ifelse(data$prin_cmp,data$rnd,0) 
        names(data)
        data<-data[order(data$rnd),]       # order data by arm
        stage1 <- glm(trt_rcv ~ rnd + x2, family=binomial, data=data)
        x.u <- residuals(stage1, type="response")
        alpha<-stage1$coeff #matrix of coeff
        v1<-vcov(stage1) #varcov matrix 
        expitWalpha<-as.vector(stage1$fitted.values)
        stage2<-glm(y ~ trt_rcv   + x2 + x.u, family=binomial, data=data)
 
        beta<-stage2$coeff #matrix of coeff
        v2<-vcov(stage2) #varcov matrix 
        expitXbeta<-as.vector(stage2$fitted.values)
        
        # Bind the  IV-variables for the rhs of the first stage GLM equation into a matrix
        # don't include the endogeneous var or xuhat
        # -- do include the IVs. **
        # -- do include a constant term **
        W<-as.matrix(cbind(data$cons,data$rnd,data$x2))
        
        # X-variables for the rhs of the second stage GLM equation into a matrix
        # don't include the endogeneous var or IV
        # -- do include xuhat **
        
        X<-as.matrix(cbind(data$x2,x.u))
        # Generate 2 matrices:
        #X0 does not include the endogeneous variable xp
        # X1 does include the endogeneous variable xp
        # Appending a constant term to the beginning of each matrix.
          X0=as.matrix(cbind(rep(1,dim(X)[1]),X)) #
        X1= as.matrix(cbind(rep(1,dim(X)[1]),data$trt_rcv,X))
        # Compute x1b1 multiplying the matrix of exogenous variables (X1) by the coefficient vectors. 
        
        x1b1=X1 %*% beta  #expitXbeta
        
  # Compute the asymptotic covariance matrix of the 2SRI estimates 
        paMu1=-beta["x.u"]*expitXbeta*expitWalpha*W[,1]
        paMu2=-beta["x.u"]*expitXbeta*expitWalpha*W[,2]
        paMu3=-beta["x.u"]*expitXbeta*expitWalpha*W[,3]

        paMu<-cbind(paMu1,paMu2,paMu3)
       
         pbMu1=inv.logit(x1b1)*X1[,1]
        pbMu2=inv.logit(x1b1)*X1[,2]
        pbMu3=inv.logit(x1b1)*X1[,3]
        pbMu4=inv.logit(x1b1)*X1[,4]
        pbMu<-cbind(pbMu1,pbMu2,pbMu3,pbMu4)
        
pbaq=t(pbMu)%*%paMu
pbbq=t(pbMu)%*%pbMu
    
D11=v1
D12=v1%*%t(pbaq)%*%solve(pbbq)
D22= solve(pbbq)%*%pbaq%*%v1%*%t(pbaq)%*%solve(pbbq)+v2
#D=D11, D12 \ t(D12), D22
D=as.matrix(rbind(as.matrix(cbind(D11,D12)),as.matrix(cbind(t(D12), D22))))
diag(D22)
 

betaiv<-stage2$coeff[2]
se.iv<-sqrt(diag(D22))[2]

        r.mu.cace[i]<-betaiv
        r.se.cace[i]<-se.iv
      } 
      
      results1<-as.data.frame(cbind(r.mu.cace,r.se.cace))
      save(results1,file=paste("binary_out/simulation_",beta.tc[trt],"/N",N[nmb],"/scenario_",scn,"/results_IV2RI_unobservedx1.Rdata",sep=""))  
      #change this for settings _wc0
      
      load(paste("binary_out/simulation_",beta.tc[trt],"/N",N[nmb],"/scenario_",scn,"/results_prin_cmp.Rdata",sep=""))  
      results1$r.cace<-results$r.cace
      result<-subset(results1,results1$r.cace<15)
      
      
      perf.meas<-performance(r.mu.cace, r.se.cace, true.delta[trt],n.sims=dim(results1)[1])
      names(perf.meas)<-c("coverage","bias","MC error bias","rmse", "median ci width") 
      write.csv(perf.meas, paste("Results/results_binout_sim",beta.tc[trt],"_N",N[nmb],"_sc",scn,"/perf_TS_fully.csv",sep=""))
      #change these for settings _wc0
      
    } 
  }
}


