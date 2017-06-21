# joint-modeling-tutorial
R code and data files associated with the joint modeling tutorial (Palestro, et al., submitted)

There are 2 R files and 2 text files corresponding to the simulation portion of the paper. The R file "r2jags_directional" contains the code used to simulate and sample from the Directed joint model in the simulation section. This file calls text file "model_directional.txt" into R. The R file "r2jags_hierarchical" contains the code used to simulate data and sample from the Hierarchical joint model in the simulation section. This file calls text file "model_hierarchical.txt" into R.
 
There are 3 data files and 1 R file needed for the application.The R file "r2jags_application" contains the code used to simulate data from the Directed joint model in the application section. The data files "Block1.Rdata" contains the behavorial data associated with the first run of the fMRI experiment, and "standardizedBOLD.Rdata" contains the neural data associated with the frist run of the fMRI (already processed). 
 
