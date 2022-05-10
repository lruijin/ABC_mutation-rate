# ABC-based estimator based on fast simulator
This folder includes the source code for ABC estimators designed for small mutation rates (< 1e-5) based on the fast simulators `fluc_exp1_rev` and `fluc_exp2_rev`. They are based on constant mutation rate assumption and piece-wise constant rate assumption, respecively. There are still three parts in the folder: I. Simulation study 1 (for the constant mutation scenario samller than 1e-5), II. Simulation study 2 (for the piece-wise constant mutation scenario), and III. Real data analysis. 

## I. Simulation study 1
* fluc_exp1_rev: simulator, generate fluctuation data for multiple parallel cultures based on constant mutation rate assumption.
* MOMMLE_fluc_exp1: calculate point estimation of MOM/MLE estimator.
* BootCI_fluc_exp1: calculate empirical confidence interval of MOM/MLE estimation using bootstrap method.
* ABC_mu1: GPS-ABC estimator based on constant mutation rate assumption.
* ABC_MCMC: ABC-MCMC estimator based on constant mutation rate assumption.
* simulation1: apply GPS-ABC estimator on a simulated dataset.

To replicate the results in the top part of Table 2, the following functions are used:
* simulation1_server: a function to run multiple replicates in simulation 1 for mutation rate from 1e-8 to 1e-4. 


## II. Simulation study 2

*Note*: The detailed model specifications are described as comments in the script.    

## III. Real data analysis
