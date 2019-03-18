# SMCMIC
Code to implement local average treatment effect via Principal Stratification mixture models  with multiple imputation estimation
#####Code for “Local average treatment effects estimation via substantive model compatible  multiple imputation”
by Karla DiazOrdaz and James Carpenter (LSHTM) 
######  
Please  address questions, comments and remarks on the code (and report bugs) to:
 	karla.diaz-ordaz@lshtm.ac.uk
as Karla DiazOrdaz is solely responsible for writing the code.
######
The code was run on R through R studio, in a MacBook Pro (R version 3.3.3 Patched (2017-08-09 r75161)) 
and also on a high performance cluster of computers (HPC)
 	
With the following packages 
* smcfcs_1.2.1
* systemfit_1.1-18
* Matrix_1.2-8
* mitools_2.3
* mice_2.30       
* boot_1.3-18     
* AER_1.2-5       
* lmtest_0.9-35    
* zoo_1.8-0        
* R2jags_0.5-7     
* rjags_4-6       
* coda_0.19-1      
* foreign_0.8-67  

########### 		
README 
###########

1. Analysis with Bayesian Mixture models
The model files “BayesMixModel_contoutcome_aux.txt” and “BayesMixModel_binoutcome_aux.txt” contain the Bayesian models corresponding to equation 2 in the manuscript, for continuous and binary outcome respectively, with the priors as specified in section 4. 
These models were run using the R2jags package from R in a HPC for parallel computing, using the scripts
“Script_BayesMixModel_Contoutcome.R” for continuous outcomes   
“Script_BayesMixModel_binoutcome_missy_auxvar.R” for binary outcomes.

The user wanting to reproduce this will need to add the model files and data files to the correct directory. A small amount of post-processing of the results needs to be done (collating the results from all 2500 simulated sets which are run in batches of 25 by the script).

The post-processing script is: "res_Bayes_bin_UX_miss_aux.R"  and "res_bayes_cont_UX_missy.R" for binary results and continuous results respectively.
These programs call the function  "perf.meas.bayes.R" to compute the performance metrics reported in the manuscript.


2. Analysis using TSLS and TSRI
The analysis of continuous outcomes using TSLS can be done using the script “”Script_2sls_contoutcome_.R“ while for binary outcome we use “Script_IV2RI_binoutcome_missoutcome.” Which calls the function performance.R  to compute the performance metrics reported in the manuscript.
A function implementing TSRI for more general use has been provided in the file “tsri.R”. The default settings are illustrated with the names corresponding to the variables in the case study presented in the paper. 

3. Analysis using SMCFCS
Sample code that uses the R package smcfcs can be Example_smcmic.R


