
library(MASS) 
library(foreign)
library(mice)
library(compiler)
source("DGP2class_contout.R")
source("DGP2class_binout.R")
DGP.bin.c <- cmpfun(DGP.bin)
DGP.cont.c <- cmpfun(DGP.cont)

n.sims=2500
#data generation
Sigma <- matrix(c(1,.3,.3,1),2,2)
set.seed= 2407
#Parameters corresponding to eq. 5 in the manuscript
#psi0=.85 for .7 proportion of compliers
#0 for .50 proportion compliers
psi1= -0.85
psi0=c(.85,.5)
#parameters in eq. 6 of the manuscript
beta.x1=-2.2
#psi2=.1 is fixed

#two sample sizes
N=c(1000, 200)
#two different  true CACE
beta.tc<-c(2,4)


##before running 
## needs to create in advance folders where the data generated will be stored 
# path of the form
#simulation_",beta.tc[trt],"/N",N[nmb],"/scenario_",scn,"/datasets/


for (trt in 1:2){
for (nmb in 1:2){
for (scn in 1:2) {
for (i in 1:n.sims){
  data <- DGP.cont(Sigma, psi0[scn], psi1, beta.tc[trt], beta.x1,  N=N[nmb], p.miss.0=0.2,theta.x2= 0.6931472)
  write.dta(data,paste("cont_out/simulation_",beta.tc[trt],"/N",N[nmb],"/scenario_",scn,"/datasets/data_",i,".dta",sep=""))
  data <- DGP.bin.c(Sigma, psi0[scn], psi1, beta.tc[trt],beta.c=0, beta.x1, N=N[nmb], p.miss.0=0.2,theta.x2= 0.6931472)
  write.dta(data,paste("binary_out/simulation_",beta.tc[trt],"/N",N[nmb],"/scenario_",scn,"/datasets/data_",i,".dta",sep=""))
#extra simulation where there is an effect of complying in controls, changing coeffeicient beta.c 
  data <- DGP.bin.c(Sigma, psi0[scn], psi1, beta.tc[trt],beta.c=beta.tc[trt]/2, beta.x1, N=N[nmb], p.miss.0=0.2,theta.x2= 0.6931472)
  write.dta(data,paste("binary_out/simulation_",beta.tc[trt],"/N",N[nmb],"/scenario_",scn,"/datasets/data_",i,"_wc0.dta",sep=""))
    }
}
}
}






