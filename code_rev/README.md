# ABC-based estimator based on fast simulator
This folder includes the source code for ABC estimators designed for small mutation rates (< 1e-5) based on the fast simulators `fluc_exp1_rev` and `fluc_exp2_rev`. They are based on constant mutation rate assumption and piece-wise constant rate assumption, respecively. There are still three parts in the folder: I. Simulation study 1 (for the constant mutation scenario samller than 1e-5), II. Simulation study 2 (for the piece-wise constant mutation scenario), and III. Real data analysis. 

## I. Simulation study 1
* fluc_exp1_rev: simulator, generate fluctuation data for multiple parallel cultures based on constant mutation rate assumption.
* MOMMLE_fluc_exp1: calculate point estimation of MOM/MLE estimator.
* BootCI_fluc_exp1: calculate empirical confidence interval of MOM/MLE estimation using bootstrap method.
* ABC_mu1: GPS-ABC estimator based on constant mutation rate assumption.
* ABC_MCMC: ABC-MCMC estimator based on constant mutation rate assumption.
* simulation1: apply GPS-ABC estimator `ABC_mu1` on a simulated dataset.

To replicate the results in the top part of Table 2, the following functions are used:
* simulation1_server: a function to run one replicate in simulation 1 for mutation rate from 1e-8 to 1e-4 using GPS-ABC estimator.
* simulation1_MCMC_server: a function to run one replicate in simulation study 1 for mutation rate from 1e-8 to 1e-4 using ABC-MCMC estimator.
* simu1_bootCI: calculates the confidence interval of MOM/MLE for simulation study 1.
* /post/simu1/*.mat : initial training samples for fitting GP surrogate model.

## II. Simulation study 2
* fluc_exp2_rev: simulator, generate fluctuation data based on piece-wise mutation rate assumption.
* ABC_mu2a: GPS-ABC estimator based on the assumption of piece-wise constant mutation rate and differential growth between mutants and normal cells.
* simulation2: apply GPS-ABC estimator defined in `ABC_mu2a` on a simulated dataset.

To replicate the results in Simulation Study 2, the following functions are used:
* simulation2_server: a function to run one replicate in simulation study 1 using GPS-ABC estimator.
* post/init_2a.mat: initial training samples for fitting GP surrogate model.
*Note*: The detailed model specifications are described as comments in the script.    

## III. Real data analysis
* ABC_mu1a: GPS-ABC estimator under the assumption of constant mutation rate and differential growth between mutants and normal cells.
* werngren_main: manually insert data (as shown in Table 5 in the manuscript), apply GPS-ABC estimators defined in ABC_mu1, ABC_mu1a and ABC_mu2a onto the datasets, and generate simulated data using estimated parameters for model selection.
* werngren_server1: a function to apply ABC_mu1 on the specified, given transition step width, initial value, expected growth rate and number of posterior samples. This function was used to test if the estimation is sensitive to the choice of initial values and hyperparameters.
* werngren_server1a: a function to apply ABC_mu1a on the specified dataset.
* werngren_server2a: a function to apply ABC_mu2a on the specified dataset.
* model_selection: summarize the numbers shown in Table 7 from the simulated data with estimated parameters.
