%%%%%%%%%%%%%%%%%% Part 0: add pathes of all functions %%%%%%%%%%%%%%%%%%%%%
cd('~/ABC/code');
addpath('~/ABC/code/kernel')
addpath('~/ABC/code/lightspeed')
addpath('~/ABC/code/optim/Matlab')

%%%%%%%%% Part 1: manually input real data and save into .mat file %%%%%%%%%
% Number of mutants of 30 samples, each has 90 cultures.
X_t = [33 18 839 47 13 126 48 80 9 71 196 66 28 17 ...
       27 37 126 33 12 44 28 67 730 168 44 50 583 23 17 24]';
% Total number of cells of 30 samples, each has 90 cultures
N_t = [1.83 1.79 1.82 1.79 2.02 2.05 1.76 1.85 2.06 2.02...
      1.9 1.9 1.9 1.9 1.9 1.9 1.9 1.9 1.9 1.9 1.9 1.9...
       1.9 1.9 1.9 1.9 1.9 1.9 1.9 1.9]'*1e8;
% Calculate the statistics as the square root of the mutant proportion
obs_X1 = sqrt(X_t(1:10)./N_t(1:10));             % of the first 10 samples
obs_X2 = sqrt(X_t(11:30)./N_t(11:30));           % of the last 20 samples
obs_X3 = sqrt(X_t./N_t);                         % of the entire 30 samples

% Calculate the moment-based estimator
MOM1 = 0.5*(1-log(mean(N_t(1:10))-mean(X_t(1:10)))/log(mean(N_t(1:10))));
MOM2 = 0.5*(1-log(mean(N_t(11:30))-mean(X_t(11:30)))/log(mean(N_t(11:30))));
MOM3 = 0.5*(1-log(mean(N_t)-mean(X_t))/log(mean(N_t)));
% Save data and initial results to data.mat
save('data.mat','obs_X1','obs_X2','obs_X3','MOM1','MOM2','MOM3');

%%%%%%%%%%%%%%%%%% Part 2: constant mutation rate estimation %%%%%%%%%%%%%%%
% model specifications:
model_spec.num_training_theta=[20,5];       % sample size for initial and refined trainings
model_spec.ksi = 0.25;                      % threshold of making wrong decision. If exceed, retrain GP model
model_spec.eps = 0;                         % discrepancy between pseudo data and real data
model_spec.M = 100;                         % samples drawn from the surrogate distribution
model_spec.N = 50000;                       % Number of posterior samples 
model_spec.a = 1;                           % Geo distribution's parameter    
model_spec.chkt = 14.56;                    % plating time, calculated as 14.56
model_spec.num_rep = 4;                     % number of replications at each training point
model_spec.time_update = 1000;              % Update the process of sampling every 1000 samples.
model_spec.L = 90;                          % Number of cultures in each sample

model_spec.init_param = [6;1e-4];           % Initial value of lengthscale and nugget
model_spec.bounds.upper = 8;                % upper bound for lengthscale parameter
model_spec.bounds.lower = 0.01;             % lower bound for lengthscale parameter

model_spec.trans_step = 5;                  % transition step width, tuned to have the acc rate between 0.2 and 0.5
model_spec.theta_old = log10(MOM3);         % Initial value of the parameter.

theta_mu = log10(MOM3);                     % prior mean
theta_sigma = 10;                           % prior sd

param_range = [-9,-4];                      % parameter range

% obtain posterior samples
POOL = parpool('local',20);
sample = ABC_mu(theta_mu, theta_sigma, param_range,obs_X3, model_spec); % obs_X3 is for the entire sample, use obs_X1 for only the first 10 samples and obs_X2 for only the last 20 samples
delete(POOL)

% calculate the acceptance rate.
acc_rate = sum(diff(sample)~=0) / model_spec.N;


%%%%%%%%%%%%%%%%%% Part 3: piecewise mutation rates estimation %%%%%%%%%%%%%%%
% reset the working space and load in data and functions again
clear;
cd('~/ABC/code');
load('data.mat');
addpath('~/ABC/code/kernel')
addpath('~/ABC/code/lightspeed')
addpath('~/ABC/code/optim/Matlab')

% model specifications
theta_mu = [log10(MOM3) log10(MOM3) 10];          % prior mean
theta_sigma = [20 20 100];                        % prior standard deviation

param_range = [-15,-2;-15,-2;0,14.56];            % range of the parameters

model_spec.num_training_theta = [1000; 10];       % initial and refinement training size
model_spec.ksi = 0.35;                            % threshold for making wrong decision, refine GP model if exceed it.
model_spec.eps = 1e-8;                            % discrepancy between pesudo data and real data
model_spec.M = 100;                               % number of samples from surrogate distribution
model_spec.N = 50000;                             % number of posterior samples
model_spec.a = 1;                                 % mean parameter in geometry distribution
model_spec.chkt = 14.56;                          % plating time
model_spec.L = 90;                                % number of cultures in one sample  
model_spec.constraint = 0;                        % If inequality assumption applies. 0 if not.

model_spec.num_rep = 3;                           % number of replications at each training point
model_spec.time_update = 500;                     % number of samples to show the process 

model_spec.init_param = [0.5 6 0.5 1]';           % Initial values of the lengthscale parameters and nugget
model_spec.bounds.upper = [8 8 2];                % upper bounds for finding optimal lengthscale parameter
model_spec.bounds.lower = [0.1 0.1 0.1];          % lower bounds for finidng optimal lengthscale parameter

model_spec.model = 2;                             % Set as 2 to use the two-stage model.
model_spec.trans_step = [2.8 0.7 1.3];            % transition step width of the three parameters
model_spec.theta_old = [log10(MOM3) log10(MOM3) 7]; % Initial values of the three parameters

% obtain posterior samples
POOL = parpool('local',32);
obs_X = sqrt(obs_X3);     % for the entire sample, the statistic we use now is the 4th root of proportion.
tStart = tic;
sample = ABC_mu2(theta_mu, theta_sigma, param_range, obs_X, model_spec);
tEnd = toc(tStart);
delete(POOL)

% You can speed up the process by calling in the initial samples provided in theta_list_full.mat.
POOL = parpool('local',32);
load('theta_list_full.mat')
obs_X = sqrt(obs_X3);
sample = ABC_mu2_init(theta_mu, theta_sigma, param_range, obs_X, theta_list, X, model_spec);
delete(POOL)
save('real_sample_xi35.mat','sample','model_spec','param_range','theta_mu');
