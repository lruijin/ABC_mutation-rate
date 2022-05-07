# Code for "Estimating mutation rates in a Markov branching process using approximate Bayesian computation"

This software package includes the source code (mostly in MATLAB) for our manuscript "Estimating mutation rates in a Markov branching process using approximate Bayesian computation". This package includes two sets of ABC estimators based on the "exact" simulator and the fast simulator, respectively. When mutation rates are at the scale less than < 10^{-5}, the exact simulators `mut_bMBP.m` and `countsizeBDtree2.m` are very slow. Estimators based on the fast simulator in folder `code_rev` are then used. There are three parts in the package: I. Simulation study 1 (for the constant mutation scenario), II. Simulation study 2 (for the piece-wise constant mutation scenario), and III. Real data analysis. 

## I. Simulation study 1
* For large mutation rate:
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
* For small mutation rate:
Codes and instructions are in the folder `code_rev`.

## II. Simulation study 2
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
