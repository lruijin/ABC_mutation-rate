# Code for "Estimating mutation rates in a Markov branching process using approximate Bayesian computation"

This software package includes the source code (mostly in MATLAB) for our manuscript "Estimating mutation rates in a Markov branching process using approximate Bayesian computation". There are three parts in the package: I. ABC-based estimator and MOM/MLE in this main folder, II Simulation study 1 (for the constant mutation scenario) and Simulation study 2 (for the piece-wise constant mutation scenario) under subfolder `simulation`, and III. Real data analysis under subfolder `real`. 

## Estimating a constant mutation rate
* ABC-based estimator with "exact" simulator:
    * **fluc_exp1.m: The "exact" simulation. Generate fluctuation data in parallel cultures from bMBP model with constant mutation rate. Used everywhere.**
    * **ABC_fluc_exp1.m: Estimate mutation rate by ABC for fluctuation data with constant mutation rate. Used in `simulation/simu1A_cascades.m` and `simulation/simu1B_cascades.m`.**
    * trainGPS.m: Train GPS model for fluctuation data with constant mutation rate. Used in ABC_fluc_exp1.m, also in demoGPS_fluc_exp1.m.
    * tnrnd.m: Generate random sample from truncated normal. Used in ABC_fluc_exp1.m.
    * mut_bMBP.m: Generate (z, x) data in a single culture from bMBP model with constant mutation rate. Used in fluc_exp1.m.
    * selectstat_fluc_exp1.m: Compare the response curves of three different summary statistics (Fig. 1).
    * demoGPS_fluc_exp1.m: Demonstrate the performance of GP regression (Fig. 2).
* ABC-based estimator with "fast" simulator
    * **fluc_exp1_rev.m: The "fast" simulator. Generate fluctuation data in parallel cultures from approximated bMBP model with constant mutation rate. Differential growth between mutants and normal cells is allowed. Used everywhere.**
    * **ABC_mu1.m: GPS-ABC estimator based on constant mutation rate assumption. Used in `simulation/simulation1_server` and `real/werngren_server1`.**
    * **ABC_MCMC.m: ABC-MCMC estimator based on constant mutation rate assumption. Used in `simulation/simulation1_MCMC_server`.**
    * **ABC_mu1a.m: GPS-ABC estimator to esimate constatnt mutation rate and relative growth rate of mutants vs. normal cells. Used in `real/werngren_server1a`.**
    * get_IniX1*.m: Generate initial training samples for fitting GP surrogate model. Used in ABC_mu1*.m. 
    * get_theta1.m and get_theta1a.m: Propose from the proposal distribution. Used in `ABC_mu1.m` and `ABC_mu1a.m`, respectively.
    * get_theta1a.m: Propose from the proposdal distributon.
    * tnorm.m: Generate random sample from truncated normal. Used in ABC_mu*.m
    * unconditionError: Calculate unconditional error of making an accpet/reject decision. Used in ABC_mu*.m
    * mut_bMBP_rev: Generate approximated (z, x) data in a single culture based on the conditional distribution of the occurring time of mutations.
* MOM/MLE estimator:
    * **MOMMLE_fluc_exp1.m: Estimate mutation rate by MOM and MLE for fluctuation data with constant mutation rate. Used in simu1B_cascades.m.**

## Estimating piece-wise constant mutation rates
* For large mutation rate `simulation2_main.m`:
    * Step 1: Generate fluctuation data from bMBP model with piecewise constant mutation rates, using function countsizeBDtree2.m.
    * Step 2: Source functions by
        * set the current folder as the working directory,
        * add the three folders in the current folder by "addpath". They contain backend functions for fitting GP model, speeding up matrix operation, and optimization.
    * Step 3: Obtain posterior samples with simulated data, using function `ABC_mu2.m`.

* For small mutation rate:
Codes and instructions are in the folder `code_rev`.

*Note*: The detailed model specifications are described as comments in the script.    

## III. Real data analysis
We apply the GPS-ABC estimator to analyze 13 datasets from a fluctuation experiment. Each dataset contains 25 parallel cultures.
Codes and instructions are in the folder `code_rev`.
