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
