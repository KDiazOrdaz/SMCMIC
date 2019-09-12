#####Code for “Local average treatment effects estimation via substantive model compatible  multiple imputation”
by Karla DiazOrdaz and James Carpenter (LSHTM) 
Biometrical Journal. 2019; 1– 15. https://doi.org/10.1002/bimj.201800345
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

Note: R has recently change the pseudo random seed generator, which may affect results. 
########### 		
README 
###########

0. Project Folder structure: In your project directory, you will need to generate the following folder structure
Two folders named "cont_out" and "binary_out". In each, the following "folder tree" should be created:
"simulation_[trt]" with [trt] taking 2 values (here 2 and 4). Inside each simulation folder,  two folders "N[nmb]" with  [nmb] taking 1000 or 200 are created. Nested within, create a further 2 folders "scenario_1" or "scenario_2" (corresponding to whether the noncompliance is 30 or 50 %). Finally, in these terminal branches, create folders "datasets" and "results".
 
At the top, in the project folder, a folder "Results" should be created, containing 16 folders, one per simulated scenario.
In the following code, these are named "results_[outc]out_sim[trt]_N[nmb]_sc[scenario]", where [outc] is either "cont" or "binary",
corresponding to continuous or binary outcomes, and the rest are as before. 
For example, one folder is "results_binout_sim2_N200_sc1".

All the scripts, functions and model files should be place directly on the Project folder and not in any subfolder.
 

1. Data generation
The simulated scenarios where created by running the file 
"script_DGP_caceBiomJ.R"
which calls the two functions "DGP2class_contout.R" to simulate continuous outcomes and "DGP2class_binout.R"  to simulate continuous outcomes.

The functions have been coded so that the options (specified in the script) correspond to the two levels of each factor variable in the simulation.

For settings where the conditional log-odds ratios is 4, some simulated datasets result in very extreme results, where the estimated LATE using the 
"true" compliance class being very large (>15). These datasets should be removed from the results in all analyses. 
We do this by running the script "princ_comp_cace.R" which saves the LATE estimated using the "true" compliance class, and saves this for all 2500 datasets in each scenario. Later those resulting in "true compliance class" LATE > 15 will be dropped. (In all cases, the results are based on at least 2000 datasets.
 
2. Analysis with Bayesian Mixture models
The model files “BayesMixModel_contoutcome_aux.txt” and "BayesMixModel_binoutcome_aux.txt” (for when there is missing data) contain the Bayesian models corresponding to equation 2 in the manuscript, for continuous and binary outcome respectively, with the priors as specified in section 4, when there is missing outcome data, and an auxiliary variable. 
These models were run using the R2jags package from R in a HPC for parallel computing, using the scripts
“Script_BayesMixModel_Contoutcome_missy.R” for continuous outcomes and  
“Script_BayesMixModel_binoutcome_missy_auxvar.R” for binary outcomes in settings with missing outcomes.

In settings with no missing data, "BayesMixModelContOutcome.txt" and “BayesMixModel_binoutcome.txt" are the Bayesian models for continuous and binary outcomes, using no covariates, which can be run using "Script_BayesMixModel_Contoutcome_fully.R" and "Script_BayesMixModel_binoutcome_fully.R", for when with continuous and binary outcome Y which is fully observed. Note: the binary scripts should be modified manually to analyse the datasets in settings where there is an effect of complying in the controls. This is easily done by modifying the datasets name (adding the _wc0 suffix), and saving the results in a different file (we also used the suffix _wc0). 

These models correspond to the (informative) priors used in the reported results (though other priors where tried for sensitivity). 

The user wanting to reproduce this will need to add the model files and data files to the correct directory.
These scripts require that the files with the Bayes models  are copied into the "terminal leaf" folder on the project, for example "CACE MI/binary_out/simulation_2/N1000/scenario_1/", which should be set as the working directory for this script manually.

When run in batches, a small amount of post-processing of the results needs to be done (collating the results from all 2500 simulated sets which are run in batches of 25 by the script).
The post-processing script is: "res_Bayes_bin.R"  and "res_bayes_cont.R" for binary and continuous results respectively.
These scripts can be modified to collate results in the no missing outcomes settings very easily (just change the name of the results files, a commented line can be uncommented to do this).

These programs call the function  "perf.meas.bayes.R" to compute the performance metrics reported in the manuscript, and save the results to the Results folder, in the corresponding subfolder (for the given scenario).

Note: the Bayesian models take longer than 3 hrs to run for the full simulation with the number of iteration specified. These could be reduced, with the potential risk of stopping the MCMC before it had reached stationarity. 
Intermediate results to reproduce the Figures in the article are given in the "Simulation_Results_Data.csv".
 
3. Analysis using TSLS and TSRI (contained in subfolder TS)
The analysis of continuous outcomes using TSLS can be done using the script “Script_TSLS_contoutcome.R”  using the package AER. 
It saves the estimates and standard errors for settings with fully observed outcomes, and then with missing outcomes, using multiple imputation by chained equations (calling the R package "mice").

The two-stage residual inclusion for binary outcome is coded manually “Script_TSRI_binoutcome_fully.R” for fully observed outcomes and Script_TSRI_binoutcome.R” for settings with missing outcome, where MI using mice is implemented. These two scripts should be modified manually to analyse the datasets in settings where there is an effect of complying in the controls. This is easily done by modifying the datasets name, and saving the results in a different file. 
The scripts contain lines that can be uncommented to do this. 

The 3 scripts call the function "performance.R" to compute the performance metrics (bias, Confidence interval coverage rate, etc) reported in the manuscript, and save these to the corresponding "Results/" subfolder.
These scripts can be done on a standard laptop (without the need to parallelise), and they are very fast. 

A function implementing TSRI for more general use has been provided in the file “tsri.R”; this implements the sandwich variance estimator of Terza 2014 (detailed in the paper). The default settings are illustrated with the names corresponding to the variables in the case study presented in the paper. 

NB. The scripts should have saved the seed for the multiple imputation, but did not. This may affect the results only a little. 
\
4. Analysis using SMCFCS
These analyses were performed in an HPC, using the scripts “smcfcs_cont.R” and “smcfcs_bin.R” for continuous and binary outcomes respectively. The user will need to run them in the same folder where the datasets folder is and NOT at the project folder root.
So, for example, change the working directory to "binary_out/simulation_",beta.tc[trt],"/N",N[nmb],"/scenario_",scn,
for example "CACE MI/binary_out/simulation_2/N1000/scenario_1/"

To analyse the scenarios for binary data where there is an effect of complying in the compliers, the analyst needs to manually modify the name of the data file in “smcfcs_bin.R” (by adding the suffix "_wc0" as before.

 
After analysing the simulated datasets, some post-processing to gather the results is needed. The scripts for this are: "res_smc_bin_miss.R" and "res_smc_bin_fully.R" (to collate results from scenarios with and without missing binary outcome), and "res_smc_cont_miss.R" and "res_smc_cont.R" for continuous outcomes. These need to be run with the root Project directory as working directory.
These programs call the function  "performance.R" to compute the performance metrics reported in the manuscript (bias, coverage rate).

To facilitate uptake, an example of a single analysis is provided in the file "Example_smcmic.R"

5. Results (tabulating and plotting)
The scripts from steps 2, 3 and 4 produce CSV files, recording the bias, coverage rate, MC error bias, RMSE and median CI width, for each method (Bayesian mixture models, TS or SMCFCS). These were collated onto a three table of results: one for continuous and two for binary outcomes, with or without an effect of compliance on the controls (collating results with the suffix "_wc0") using the script "results_tables.R".

Since both the Bayesian and the SMCFCS analysis take considerable time when running in all 2500 replicates for all the scenarios, we include here tables of intermediate results to allow to reproduce the plots.

The plots were produce with the code contained in "Results_plots.R"
The intermediate table of results, used for the Figures reported in the published paper are "cont_results_intermediate.csv" for the continuous outcomes and "bin_results_intermediate.csv" and "bin_results_wc0_intermediate.csv" for the binary outcome results (with no effect of complying in the control arm and with a positive effect, respectively). 
These tables of results should be saved in the "Results" folder of the Project, and used with "Results_plots.R" to produce Figure 1, 2 and 3 respectively on the manuscript. 

 


