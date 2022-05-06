%%%%%%%%%%%%%%%%%%%%%%%%% Part 1: Generate data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%% Step 1: Set the hyperparameters %%%%%%%%%%%%%%%%%%%
  % Exponential rate parameter a, typically 1 
    a = 1; % may change a to other values
    p = 1e-4;
    c = 100;
    Z0 = 1; % may change Z0 to other values
    myfun = @(t, Z0, a, p, c) Z0 * (exp(a * t) - exp(a * t * (1 - 2 * p))) - c;
    fun = @(t) myfun(t, Z0, a, p, c);
    tp = fzero(fun, 20);
    L = 1;
    nmcs=[10 50 100];
%%%%%%%%%%%%%%%%%%%%%%%%% Step 2: Generate fluctuation data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Get multiple local workers
POOL = parpool('local',10);
nSimu = 100;
% Set seeds
for seed = 1:nSimu;
  if mod(seed, 5) == 0
        fprintf(['seed=', int2str(seed),'\n']); % print  the process
      end
  for k = 1:3
    nmc = nmcs(k);
    startTime = tic;
    s = RandStream('mt19937ar','Seed',seed);
    RandStream.setGlobalStream(s);
    % Initialize number of total cells and number of mutants
    Nt = NaN(nmc,1); 
    Xt = NaN(nmc,1);
    % Initialize the vector for saving observed data: the statistics
    obs_X = NaN(nmc,1);

    parfor j = 1 : nmc
      [tmp1, tmp2] = mut_bMBP_rev(Z0, a, 1, p, tp);
      Nt(j,1) = tmp1;%%
      Xt(j,1) = tmp2;
      obs_X(j,1) = sqrt(tmp2/tmp1);%% Statistic is the 4th root of mutant proportion
    end
    % calculate the moment-based estimator
    [MOM MLE] = MOMMLE_fluc_exp1(Nt, Xt);
    
    % save simulated data, including Nt, Xt, obs_X and MOM
    duration = toc(startTime);
      save(strcat('/data/ZChenLab/Lu/ABC/data/simu1/p4/Y',num2str(k),'_', ...
      num2str(seed),'.mat'),'Nt','Xt','obs_X','MOM','MLE','tp','duration');
  end
end
delete(POOL);

%%%%%%%%%%%%%%%%%% Part 2: Estimating piecewise constant mutation rates %%%%%%%%%%%%%%%%
  % including folders containing functions to be used
% The main folder contains the functions to be directly used
clear;
cd('~/ABC/code_rev')
% This folder contains the functions needed for fitting GP models
addpath('~/ABC/code/kernel')
% This folder contains functions for quick matrix operations
addpath('~/ABC/code/lightspeed')
% This folder contains functions for optimization, which is for estimating hyperparameters in GP models
addpath('~/ABC/code/optim/Matlab')

%%%%%%%%% Step 1: Model specifications for the ABC estimator %%%%%%%%%%%%%%%%%%%
theta_sigma = 20;        % sd of prior distribution
p_range = [-10 -6];

% model specifications
model_spec.num_training_theta = [1500; 10];    % Number of settings to train and refine GP model 
model_spec.ksi = 0.25;                         % The threshold of incorrect decision to refine model
model_spec.eps = 0;                         % Discrepancy allowed between the true data and simulated data
model_spec.M = 100;                            % Number of samples drawn from the surrogate model
model_spec.N = 3000;                          % Total number of samples to be drawn from the posterior
model_spec.a = 1;                              % Hyperparameter of the geo dist, known, same as the one to generate data
model_spec.Z0 = 1;                              % Number of cultures in one sample, 1 culture per sample in simulation.
model_spec.saveSample = 0;
model_spec.num_rep = 10;                       % The number of replications at each setting 
model_spec.time_update = 500;                 % The number to update process onto the screen
model_spec.constraint = 0;
%kern_param = [0.01,0.4];
    model_spec.init_param = [1 1]';            % Hyperparameter for GP model, in the order of lengthscale parameters and the nugget parameter 
    model_spec.bounds.upper = 4;             % Upper bounds for the three lengthscale parameters
    model_spec.bounds.lower = 0.1;       % Lower bounds for the three lengthscale parameters
    % case 1
    %model_spec.trans_step = [0.8 3.5 3.5 4];         % The transition steps of the MCMC algorithm 
    %model_spec.theta_old = [log10(0.8) -9 -8 14];    % t_d = 4  % The initial value of the three parameters
    % case 2
    model_spec.trans_step = 0.35;         % The transition steps of the MCMC algorithm 
    model_spec.theta_old = -8;    % t_d = 4  % The initial value of the three parameters
    
    model_spec.init = 'post/simu1/init1_p8.mat';
%%%%%%%%% Step 2: Piece-wise constant mutation rates estimation %%%%%%%%%%%%%%%%%%%

%POOL = parpool('local',20);                  % always keep this as ABC_mu* functions contain parallel computations for fitting GP model.
for seed = 2:100                             % change seed range to run less cases, each takes around 40 minutes.
    load(strcat('/data/ZChenLab/Lu/ABC/data/simu1/p8/Y1_',num2str(seed),'.mat'))
    model_spec.chkt = tp;
    theta_mu = log10(MOM);   % Initial value corresponding to MOM estimator`
    [sample, runningTime, initTime] = ABC_mu1(theta_mu, theta_sigma, 1, p_range, sqrt(obs_X), model_spec);
    acc_rate = sum(diff(sample(:,1))~=0) / model_spec.N;
    % save the sample with its model specification, acceptance rate and running time.
    save(strcat('/data/ZChenLab/Lu/ABC/result/post/sample1_',num2str(seed),'.mat'), 'model_spec','theta_mu','theta_sigma','sample','runningTime','dt_range','p_range');
end
delete(POOL);
