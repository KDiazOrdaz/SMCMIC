rm(list=ls())

# This runs the models where there is no missing data
#need to change the directory to where the dataset folder is
#"cont_out/simulation_",beta.tc[trt],"/N",N[nmb],"/scenario_",scn,
# for example 
## "cont_out/simulation_2/N1000/scenario_1/"


# get running counter from the system
taskID <- as.numeric(Sys.getenv("SGE_TASK_ID"))
print(taskID) 

# set specific values for the run 
everyout <- 25 # export results every 'everyout' runs 

require(foreign)
library(R2jags)
library(MASS) 


n.sims <- 2500
count <- 0


###################################################################################

r.mu.cace <- r.se.cace <-r.lb.cace<-r.ub.cace<-r.med.cace<-numeric(n.sims)
# Parameters: pick the parameters to monitor
params <- c("a")
thin <- 1
# Model: specify the file with the model
mfile <- "../../../../BayesMixModelContOutcome.txt"




while (!file.exists(paste("results/BayesUnobsX.fully.",taskID,"to",(taskID+everyout-1),".RData",sep="")))
{
  for (i in taskID:(taskID+everyout-1)){
    count<-count+1
    cat("loop id= ",i,"\n")
    
    dta.df<-read.dta(paste("datasets/data_",i,".dta",sep=""))
    
    
    N.obs <- dim(dta.df)[1]
    y <- dta.df$y
    rnd <- dta.df$rnd
    cmp <- dta.df$cmp
    
    
    data.list <- list("N.obs", "y",  "rnd", "cmp")
    
    ##
    
    ###################################################################################
    
    
    
    ###########################################################################
    ###  Initial run to start chains and find out possible autocorrelation  ###
    ###########################################################################
    
    set.seed(100928)
    jini <- jags(data=data.list, inits=NULL, params, DIC=FALSE, n.chains=2, 
                 n.iter=10000, n.burnin=5000, n.thin=1, model.file=mfile)
    # Check autocorrelation
    jmcmc <- as.mcmc(jini)
    
    thin <- 1
    iter <- thin * 10000
    # Initialise counters
    jupd <- jini
    jupd <- update(jupd, n.iter=iter, n.thin=thin)
    jmcmc <- as.mcmc(jupd)
    # Summary statistics
    
    # saving the Summary statistics
    trt.effect<-t(summary(jmcmc)$statistics)[,"a[3]"]
    mean.cace<-trt.effect[1]
    se.cace<-trt.effect[2]
    quant.cace<-t(summary(jmcmc)$quantiles)[,"a[3]"]
    lb.cace<-quant.cace[1]
    med.cace<-quant.cace[3]
    ub.cace<-quant.cace[5]
    
    
    
    r.mu.cace[[count]]<-mean.cace
    r.se.cace[[count]]<-se.cace
    r.lb.cace[[count]]<-lb.cace
    r.ub.cace[[count]]<-ub.cace
    r.med.cace[[count]]<-med.cace
    
  } 
  
  
  res<-as.data.frame(cbind(r.mu.cace,r.se.cace,r.lb.cace,r.ub.cace,r.med.cace))
  save(res,file=paste("results/BayesUnobsX.fully.",taskID,"to",(taskID+everyout-1),".RData",sep=""))
  
  
  
}

