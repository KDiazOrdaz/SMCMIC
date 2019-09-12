DGP.cont<-function(
  sigma.mat= matrix(c(1,.3,.3,1),2,2),
  psi0=.85, #psi0 is the basic logodds of being a complier, 0 for 50% prop complier
  psi1=-2 ,   #psi1 is the coef of x1 on prin.cmp psi2 is for x2
  beta.tc=2, # Association between treatment and outcome in compliers
  beta.x1=-1, # Association between indiv var x1 and outcome
  beta.x2=0.5, # Association between indiv var x2 and outcome
  N=1000,
  p.miss.0=0.2, # % of expected missing outcome 
  #p.miss.x=0.1,	# % of expected missing outcome	
  theta.x2	 = 0.6931472,		# Association between indiv var and p.miss by enpoint
  theta.trt = 0		# Association between treatment and p.miss by endpoint
  )
  {
X<-mvrnorm(n=N, mu=c(0,0),Sigma=sigma.mat)
x1<-c(X[,1])
x2<-c(X[,2])
cons<-rep(1,N)
rnd<-rbinom(N,1,.5) # randomised allocation vector
eta.cmp<-psi0 + psi1*x1 #compliance may depend on covariates, in this case only x1
p.cmp <- 1/(1 + exp(-(eta.cmp)))
prin.cmp <- rbinom(N, 1,p.cmp)
#observed compliance
cmp<-ifelse(rnd>0,prin.cmp,NA)
##outcome beta0=0, no effect of complying in the control group, no effec of randomisation
eta.y<-beta.x1*x1+ beta.x2*x2+ beta.tc*prin.cmp*rnd +.75*prin.cmp   
y<-rnorm(N, eta.y,1)
#introduce MAR missing outcomes
theta.miss.0 <-log(p.miss.0/(1-p.miss.0))
theta.miss.y<- theta.miss.0 +  theta.x2*x2	+ theta.trt*rnd #but no interaction between x2 and rnd
  p.miss.y <- exp(theta.miss.y)/(1+exp(theta.miss.y))

  miss.y<-rbinom(N,1,prob=p.miss.y) #generate  missing outcome indicator
y.miss<-ifelse(miss.y==0,y,NA)
cons<-rep(1,N)
data<-as.data.frame(cbind(y, y.miss, x1, x2, rnd, prin.cmp, cmp, cons))
return(data)
}


