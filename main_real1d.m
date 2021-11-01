clear;
workpath = 'C:\Users\lruijin\Documents\MATLAB\ABC_wu\Code';
resultpath = 'C:\Users\lruijin\Documents\MATLAB\ABC_wu\Result';
addpath(workpath);
addpath('C:\Users\lruijin\Documents\MATLAB\ABC_wu\Code\kernel');

% X_t = [33 18 839 47 13 126 48 80 9 71 196 66 28 17 ...
%        27 37 126 33 12 44 28 67 730 168 44 50 583 23 17 24]';
% N_t = [1.83 1.79 1.82 1.79 2.02 2.05 1.76 1.85 2.06 2.02...
%        1.9 1.9 1.9 1.9 1.9 1.9 1.9 1.9 1.9 1.9 1.9 1.9...
%        1.9 1.9 1.9 1.9 1.9 1.9 1.9 1.9]'*1e8;
% obs_X1 = sqrt(X_t(1:10)./N_t(1:10));
% obs_X2 = sqrt(X_t(11:30)./N_t(11:30));
% obs_X3 = sqrt(X_t./N_t);
% 
% MOM1 = 0.5*(1-log(mean(N_t(1:10))-mean(X_t(1:10)))/log(mean(N_t(1:10))));
% MOM2 = 0.5*(1-log(mean(N_t(11:30))-mean(X_t(11:30)))/log(mean(N_t(11:30))));
% MOM3 = 0.5*(1-log(mean(N_t)-mean(X_t))/log(mean(N_t)));
% save('data.mat','obs_X1','obs_X2','obs_X3','MOM1','MOM2','MOM3');

model_spec.num_training_theta=[20,5];
model_spec.ksi = 0.25;
model_spec.eps = 0;
model_spec.M = 100;
model_spec.N = 10000;
model_spec.a = 1;
model_spec.chkt = 19;
model_spec.num_rep = 4;
model_spec.time_update = 1000;
model_spec.L = 90;

model_spec.init_param = [6;1e-4];
model_spec.bounds.upper = 8;
model_spec.bounds.lower = 0.01;

%trans_step = [0.08,0.08,0.08];
model_spec.trans_step = 14;
model_spec.theta_old = log(MOM3);

theta_mu = log(MOM3);
theta_sigma = 10;

param_range = [-20,-10];
POOL = parpool('local',20);
sample = ABC_mu(theta_mu, theta_sigma, param_range,obs_X3, model_spec);
delete(POOL)
acc_rate = sum(diff(sample)~=0) / model_spec.N;

%% The Generalized birth-death process model for estimating real data
theta_mu = [log(MOM3) log(MOM3) 2.3];
theta_sigma = [20 20 100];

param_range = [-20,-10;-20,-10;-0.1,2.95];
model_spec.num_training_theta = [150; 10];
model_spec.ksi = 0.3;
model_spec.eps = 1e-8;
model_spec.M = 100;
model_spec.N = 50000;
model_spec.a = 1;
model_spec.chkt = 19;

model_spec.num_rep = 3;
model_spec.time_update = 5000;

%kern_param = [0.01,0.4];
model_spec.init_param = [6 6 1 1]';
model_spec.bounds.upper = [8 8 2];
model_spec.bounds.lower = [0.1 0.1 0.1];

model_spec.model = 2;
%trans_step = [0.08,0.08,0.08];
model_spec.trans_step = [11 11 1];
model_spec.theta_old = [log(MOM3) log(MOM3) 2.3];

POOL = parpool('local',20);
ABC_mu2(theta_mu, theta_sigma, param_range,obs_X3, model_spec);
delete(POOL)