##### imputing latent compliance classes ####
#### by using SMCFCS package 
#### and computing CACE estimates by using max likelihood 

library(mice)
library("mitools")
library(smcfcs)
numit=250
rjlimit=5000
m<-100

### y outcome
### rnd is the randomised treatment indicator
### cmp is the compliance class, partially latent in the data

data<-read.dta(paste("datasets/data_",i,".dta",sep=""))
originaldata<-as.data.frame(cbind(data$y,data$rnd,data$cmp))
### define the interaction term, whose coefficient will identify the cace
originaldata$cace<-data$rnd*data$cmp
names(originaldata)<-c("y", "rnd","cmp","cace" )
### no missing continuous outcomes ###
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
    
##### binary outcomes , fully observed ###
    ## only change the type of substantive model type smtype
imp<-smcfcs(originaldata,smtype = "logistic",
      smformula="y~  cmp + cace +",
      method=c("","","", "","logreg","cmp*rnd"),
      predictorMatrix=NULL,
      m=m,
      numit=numit,
      rjlimit=rjlimit,
      noisy=FALSE) 

impobj <- imputationList(imp$impDatasets)
models <- with(impobj, glm(y~  cmp + cace, family=binomial()))
mean.cace<-summary(MIcombine(models))[3,1]
se.cace<-summary(MIcombine(models))[3,2]

    
###missing continuous outcomes ###
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

### binary missing outcomes #####
imp.miss<-smcfcs(originaldata.miss,smtype = "logistic",
                 smformula="y~  x2 + cmp + cace",
                 method=c("","","", "","logreg","cmp*rnd"),
                 predictorMatrix=NULL,
                 m=m,
                 numit=numit,
                 rjlimit=rjlimit,
                 noisy=FALSE) 


impobj.miss <- imputationList(imp.miss$impDatasets)
models.miss <- with(impobj.miss, glm(y~  cmp + cace, family=binomial()))

mean.cace.miss<-summary(MIcombine(models.miss))[3,1]
se.cace.miss<-summary(MIcombine(models.miss))[3,2]

    
