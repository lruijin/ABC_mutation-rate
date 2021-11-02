# Estimate mutation rates in a Markov branching process using approximate Bayesian computation

Through modeling fluctuation experimental data by a two-type Markov branching process, we use approximate Bayesican computation (ABC) to estimate parametres for time-varying mutation rates. Though complicated dynamic structure could be easily incorporated, codes in this repository are for constant and piece-wise constant mutation rate functions only.

## Simulation study 2
``simulation2_main.m`` is the main script to reproduce the results in the second simulation study. It contains the following steps:
  1. Generate data from two-type Markov branching process with piecewise constant mutation rates, using function ``countsizeBDtree2``.
  2. Add functions in by
      - set the current folder as the working directory
      - add the three folders in the current folder by ``addpath``. They contain backend functions for fitting GP model, speeding up matrix operation, and optimization.
  3. Obtain posterior samples with simulated data, using function ``ABC_mu2``.

*Note*: 
  - The detailed model specifications are described as comments in the script.    
  - The current file generates 100 datasets (seed from 1 to 100) but only conducts estimation for the first dataset. The upper bound of the for loop could be changed to play around with multiple datasets.

## Real data analysis
``real_main.m`` is the main script to reproduce results of the real data analysis. It contains the following steps:
1. Manually input the benchmark dataset from a bacterioa fluctuation experiment. It contains 30 parallel cultures. Each culture was incubated from 90 bacterium cells. Thus, we:
    - save the total number of cells of each culture in 'N_t' and the number of mutants in 'X_t'. 
    - calculate summary statistics for the first 10 cultures, the last 20 cultures and the entire 30 cultures and save them into 'obs_X1', 'obs_X2' and 'obs_X3' respectively.
    - calculate the moment-based estimators for the first 10 cultures, the last 20 cultures and the entire 30 cultures and save them into 'MOM1', 'MOM2' and 'MOM3' respectively.
2. Estimate mutation rate under the constant mutation assumption, using function `ABC_mu`. Detailed model specifications are given as comments in the script.
3. Estimate mutation rate under the piece-wise constant mutations assumption, using function `ABC_mu2`. Detailed model specifications are given as comments in the script.

*Note*:
  - Preparing initial samples is time-consuming. Thus, we provide the initial training set in `theta_list_full.mat`. Function `ABC_mu2_init` is then used to obtain the posterior samples, where preparing initial samples for GP training is skipped. 
