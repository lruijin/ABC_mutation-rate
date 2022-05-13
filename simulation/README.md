# Codes for simulation study 1 and simulation study 2
This folder includes the source code for the two simulation studies

## I. Simulation study 1
* Mutation rate larger than 1e-4 (etimators with "exact" simulator)
    * simu1B_cascades.m: Simulation study 1 (for GPS-ABC). *Note*: due to intensive computations, simu1B_cascades.m and simu1A_cascades.m run on servers. 
    * simu1A_cascades.m: Simulation study 1 (for ABC-MCMC). *Warning*: running ABC-MCMC in this simulation study takes significant amount of time.
* Mutation rate smaller than 1e-4 (estimators with "fast" simulator)
    * simulation1: Generate 100 simulated dataset using the "fast" simulator and apply GPS-ABC estimator `ABC_mu1` on a simulated dataset on one of them.
    * simulation1_server: a function to run one replicate in simulation 1 for mutation rate from 1e-8 to 1e-4 using GPS-ABC estimator.
    * simulation1_MCMC_server: a function to run one replicate in simulation study 1 for mutation rate from 1e-8 to 1e-4 using ABC-MCMC estimator.
* MOM/MLE estimator's confidence interval:
    * BootCI_fluc_exp1: calculate empirical confidence interval of MOM/MLE estimation using bootstrap method.
    * simu1_bootCI: calculates the confidence interval of MOM/MLE for simulation study 1.
    * ../traning/simu1/*.mat : initial training samples for fitting GP surrogate model.

*Note*: To run 100 replicates on server, detailed instructions are given in the header comments of the scripts.

## II. Simulation study 2
* simulation2: Generate 100 simulated datasets using the "fast" simulatorapply GPS-ABC estimator defined in `ABC_mu2a` on a simulated dataset.
* simulation2_server: a function to run one replicate in simulation study 1 using GPS-ABC estimator.
* ../training/init_2a.mat: initial training samples for fitting GP surrogate model.

*Note*: To run 100 replicates on server, detailed instructions are given in the header comments of the scripts.
