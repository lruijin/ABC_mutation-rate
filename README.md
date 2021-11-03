# Code for "Estimating mutation rates in a Markov branching process using approximate Bayesian computation"

This software package includes the source code (mostly in MATLAB) for our manuscript "Estimating mutation rates in a Markov branching process using approximate Bayesian computation". There are three parts in the package: I. Simulation study 1 (for the constant mutation scenario), II. Simulation study 2 (for the piece-wise constant mutation scenario), and III. Real data analysis.

## I. Simulation study 1
* simu1B_cascades.m: Simulation study 1 (for GPS-ABC). *Note*: due to intensive computations, simu1B_cascades.m and simu1A_cascades.m run on servers.
* simu1A_cascades.m: Simulation study 1 (for ABC-MCMC). *Warning*: running ABC-MCMC in this simulation study takes significant amount of time.
* MOMMLE_fluc_exp1.m: Estimate mutation rate by MOM and MLE for fluctuation data with constant mutation rate. Used in simu1B_cascades.m.
* ABC_fluc_exp1.m: Estimate mutation rate by ABC for fluctuation data with constant mutation rate. Used in simu1A_cascades.m and simu1B_cascades.m.
* trainGPS.m: Train GPS model for fluctuation data with constant mutation rate. Used in ABC_fluc_exp1.m, also in demoGPS_fluc_exp1.m.
* tnrnd.m: Generate random sample from truncated normal. Used in ABC_fluc_exp1.m.
* fluc_exp1.m: Generate fluctuation data in parallel cultures from bMBP model with constant mutation rate. Used everywhere.
* mut_bMBP.m: Generate (z, x) data in a single culture from bMBP model with constant mutation rate. Used in fluc_exp1.m.
* selectstat_fluc_exp1.m: Compare the response curves of three different summary statistics (Fig. 1).
* demoGPS_fluc_exp1.m: Demonstrate the performance of GP regression (Fig. 2).

## II. Simulation study 2
* simulation2_main.m: Simulation study 2 (for GPS-ABC only). It contains the following steps:
1. Generate fluctuation data from bMBP model with piecewise constant mutation rates, using function countsizeBDtree2.m.
2. Source functions by
    * set the current folder as the working directory
    * add the three folders in the current folder by "addpath". They contain backend functions for fitting GP model, speeding up matrix operation, and optimization.
3. Obtain posterior samples with simulated data, using function ABC_mu2.m.
*Note*: The detailed model specifications are described as comments in the script.    

## III. Real data analysis
* real_main.m: Real data analysis (for GPS-ABC only). It contains the following steps:
1. Manually input the benchmark dataset from a bacterioa fluctuation experiment. It contains 30 parallel cultures. Each culture was incubated from 90 bacterium cells. Thus, we:
    - save the total number of cells of each culture in 'N_t' and the number of mutants in 'X_t'. 
    - calculate summary statistics for the first 10 cultures, the last 20 cultures and the entire 30 cultures and save them into 'obs_X1', 'obs_X2' and 'obs_X3' respectively.
    - calculate the moment-based estimators for the first 10 cultures, the last 20 cultures and the entire 30 cultures and save them into 'MOM1', 'MOM2' and 'MOM3' respectively.
2. Estimate mutation rate under the constant mutation assumption, using function ABC_mu.m. Detailed model specifications are given as comments in the script.
3. Estimate mutation rates under the piece-wise constant mutations assumption, using function ABC_mu2.m. Detailed model specifications are given as comments in the script.

*Note*: Preparing initial samples is time-consuming. Thus, we provide the initial training set in theta_list_full.mat. Function ABC_mu2_init.m is then used to obtain the posterior samples, where preparing initial samples for GP training is skipped. 
