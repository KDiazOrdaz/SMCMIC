### Function for 2SRI with  SE ###
# requires library(mice) 
tsri <- function(data, #the multiply imputed datasets stacked and with identifier mj
                 A = c("trt_rcv"), #name of endogeneous exposure
                 instrument = c("treat"), #name of instrument
                 y_form  = paste("Y ~ trt_rcv + hads_dep_base +  centre + Gender + age + employed_ftedu ", collapse=""), #substantive model
                 Y_covnames = c("hads_dep12m"), #name of outcome in dataset
                 M = 50 #number of imputations
                 )
{
#select the numbe column in the dataframe that contains the endogeneous regressor
Acol <- (1:ncol(data))[colnames(data) %in% A]
Zcol <- (1:ncol(data))[colnames(data) %in% instrument]

exog.vars<-setdiff(attr(terms(as.formula(y_form)), "term.labels"),  A)
#exxog variables needed to specify the implied first stage
a_form  = paste(names(data)[Acol], "~",paste(c(instrument,exog.vars),collapse=" + "), collapse = "") 

#y_form doesn't contain residuals, need to add
stage2_form <- paste(y_form, "x.u",sep=" + ")

#adds 2 extra columns one for intercept, the other for residual from 1st stage inclusion
coef.iv<-var.coef<-matrix(NA, nrow=M, ncol=2+length(attr(terms(as.formula(y_form)), "term.labels")) )



for (j in 1:M){
  data<-subset(data_mi, data_mi$mj==j)  
  cons<-data$cons<-rep(1,dim(data)[1])
  #include the Y column
  Ycol <- (1:ncol(data))[colnames(data) %in% Y_covnames]
  data$Y<- data[,Ycol]
  #Y<- data[,Ycol]
  
  
  stage1 <- glm(as.formula(a_form) , family=binomial, data=data)
  data$x.u<- x.u <- residuals(stage1, type="response")
  
  #mata: alpha=st_matrix("e(b)")'
  #        mata: v1=st_matrix("e(V)")
  alpha<-stage1$coeff #matrix of coeff
  v1<-vcov(stage1) #varcov matrix 
  expitWalpha<-as.vector(stage1$fitted.values)
  #stage1$linear.predictors
  
  #Apply GLM for the 2SRI second stage. 
  
  stage2<-glm(as.formula(stage2_form), family=binomial, data=data)
  
  beta<-stage2$coeff #matrix of coeff
  v2<-vcov(stage2) #varcov matrix 
  expitXbeta<-as.vector(stage2$fitted.values)
  #length(beta)
  # Bind the  IV-variables for the rhs of the first stage GLM equation into a matrix
  # don't include the endogeneous var or xuhat
  # -- do include the IVs. **
  # -- do include a constant term **
  Exo.covcols <- (1:ncol(data))[colnames(data) %in% exog.vars]
  Exo.matrix <- data.frame(data[, Exo.covcols])
  
  W<-data.matrix(cbind(cons,data[,Zcol],data[, Exo.covcols]))
  #dim(W)
  # X-variables for the rhs of the second stage GLM equation into a matrix
  # don't include the endogeneous var or IV
  # -- do include xuhat **
  
  X<-data.matrix(cbind(Exo.matrix, x.u))
  # Generate 2 matrices:
  #X0 does not include the endogeneous variable xp
  # X1 does include the endogeneous variable xp
  # Appending a constant term to the beginning of each matrix.
  X0=as.matrix(cbind(cons,X)) #
  X1= as.matrix(cbind(cons,data$trt_rcv,X))
  
  
  # Compute x1b1 multiplying the matrix of exogenous variables (X1) by the coefficient vectors. 
  
  x1b1=X1 %*% beta  #expitXbeta
  
  # Compute the asymptotic covariance matrix of the 2SRI estimates 
  #
  #need to find the correct way to multiply
  #beta["x.u"]*expitXbeta*expitWalpha*W
  paMu <- matrix(data=rep(NA,dim(W)[1]*dim(W)[2]), ncol=dim(W)[2])
  for (k in 1:dim(W)[2]){
    paMu[,k]=-beta["x.u"]*expitXbeta*expitWalpha*W[,k]}
  pbMu <- matrix(data=rep(NA,dim(X1)[1]*dim(X1)[2]), ncol=dim(X1)[2])
  
  for (k in 1:dim(X1)[2]){
    pbMu[,k]=plogis(x1b1)*X1[,k]}
  
  
  pbaq=t(pbMu)%*%paMu
  pbbq=t(pbMu)%*%pbMu
  
  D11=v1
  D12=v1%*%t(pbaq)%*%solve(pbbq)
  D22= solve(pbbq)%*%pbaq%*%v1%*%t(pbaq)%*%solve(pbbq)+v2
  #D=D11, D12 \ t(D12), D22
  D=as.matrix(rbind(as.matrix(cbind(D11,D12)),as.matrix(cbind(t(D12), D22))))
  diag(D22)
  
  
  
  coef.iv[j,]<-stage2$coeff
  colnames(coef.iv) <- names(stage2$coeff)
  var.coef[j,]<-diag(D22)
  colnames(var.coef) <- names(diag(D22))
  
  
}



###### Apply Rubin's rules if M >1 ####
if (M>1){
mean.cace <-apply(coef.iv, 2, mean)
var.cace <-apply(var.coef,2,mean)
with.var.cace<- var.cace   # within-imputation variance
bet.var.cace<-(apply(coef.iv, 2, sd))^2

se.cace<- sqrt(with.var.cace+(1+1/M)*bet.var.cace)
#normal based confidence intervals 
lb.cace<- mean.cace -1.96 * se.cace
ub.cace<- mean.cace +1.96 * se.cace
}
else{
  mean.cace <-stage2$coeff
  se.cace<- sqrt(diag(D22))
  lb.cace<- mean.cace -1.96 * se.cace
  ub.cace<- mean.cace +1.96 * se.cace
  }
return(list(est.cace=mean.cace,se.cace=se.cace, lb.cace=lb.cace, ub.cace=ub.cace))

}
