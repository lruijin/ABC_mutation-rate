clear;
workpath = 'C:\Users\lruijin\Documents\MATLAB\ABC_wu\Code';
resultpath = 'C:\Users\lruijin\Documents\MATLAB\ABC_wu\Result';
addpath(workpath);
addpath('C:\Users\lruijin\Documents\MATLAB\ABC_wu\Code\kernel');
nmc = 1; mu = exp(-18); chkt = 19; a = 1;

seed = 2;
s = RandStream('mt19937ar','Seed',seed);
RandStream.setGlobalStream(s);

obs_X = NaN(nmc,1); 
Xt = NaN(nmc,1);
Nt = NaN(nmc,1);
tic;
for j = 1 : nmc
    [temp1, temp2] = countsizeBDtree(a, mu, chkt);%%
    obs_X(j, 1) = sqrt(temp2/temp1);%%
    Xt(j) = temp2;
    Nt(j) = temp1;
end
toc;
MOM = 0.5*(1-log(mean(Nt)-mean(Xt))/log(mean(Nt)));
MLE_eqn = @(x) (1-2*x)*(mean(Nt)-mean(Xt))+x-(1-x)*(mean(Nt))^(1-2*x);
MLE = log(fsolve(MLE_eqn,MOM));

save('Y.mat','obs_X');

model_spec.num_training_theta=[20,5];
model_spec.ksi = 0.25;
model_spec.eps = 0;
model_spec.M = 100;
model_spec.N = 5000;
model_spec.a = a;
model_spec.chkt = chkt;
model_spec.num_rep = 2;
model_spec.time_update = 2;

model_spec.init_param = [6;1e-4];
model_spec.bounds.upper = 8;
model_spec.bounds.lower = 0.01;

%trans_step = [0.08,0.08,0.08];
model_spec.trans_step = 1;
model_spec.theta_old = log(MOM);

theta_mu = log(MOM);
theta_sigma = 100;

param_range = [-20,-10];
sample = ABC_mu(theta_mu, theta_sigma, param_range,obs_X, model_spec);
acc_rate = sum(diff(sample)~=0) / model_spec.N;
%save(strcat(resultpath,'\sample.mat'),'model_spec','obs_X','sample','acc_rate');

% acceptance rate
sum(diff(sample(5000:10000))~=0) / 5000 % acceptance rate: 0.3418
mean(sample(5000:10000))
% plot the trace plot
figure;
plot(1:10000,sample)
xlim([0,15000]);
hold on 
plot(1:15000,repmat(-6.5,15000,1));
% plot histogram after burn-in
figure
histogram(sample(5000:15000),'Normalization','probability')
hold on;
plot([-6.5,-6.5],[0,0.05])


%% test the 2nd case - for three unkown parameters
% generate the data and plot
a = 1; mu_vec = exp([-8 -4]); t_d = 4; chkt = 9; nmc = 100;

seed = 1;
s = RandStream('mt19937ar','Seed',seed);
    RandStream.setGlobalStream(s);

Nt = NaN(nmc,1); 
Xt = NaN(nmc,1);
obs_X = NaN(nmc,1);

tic;
for j = 1 : nmc
    [temp1, temp2] = countsizeBDtree2(a, mu_vec,t_d, chkt);%%
    Nt(j,1) = temp1;%%
    Xt(j,1) = temp2;
    obs_X(j,1) = sqrt(sqrt(temp2./temp1));%%
end

MOM = 0.5*(1-log(mean(Nt)-mean(Xt))/log(mean(Nt)));
theta_mu = [-6 -2 1.5];
theta_sigma = [20 20 100];

param_range = [-9,-2;-9,-1;-0.1,2.2];

model_spec.num_training_theta = [150; 10];
model_spec.ksi = 0.3;
model_spec.eps = 1e-8;
model_spec.M = 100;
model_spec.N = 50000;
model_spec.a = 1;
model_spec.chkt = 9;

model_spec.num_rep = 10;
model_spec.time_update = 1000;

%kern_param = [0.01,0.4];
model_spec.init_param = [6 6 1 1];
model_spec.bounds.upper = [8 8 2];
model_spec.bounds.lower = [0.1 0.1 0.1];

model_spec.model = 2;
%trans_step = [0.08,0.08,0.08];
model_spec.trans_step = [6 1 0.5];
model_spec.theta_old = [-8 -4 1.3];

sample_t4 = ABC_mu2(theta_mu, theta_sigma, param_range,obs_X, model_spec);
acc_rate = sum(diff(sample)~=0) / model_spec.N;
save('../Result/sample_3d.mat','sample','model_spec','obs_X','a','chkt','mu_vec','t_d','theta_mu','theta_sigma')
load('../Result/sample_3d.mat')
figure
subplot(1,2,1)
histogram(exp(sample(10000:30000,1)),'Normalization','probability')
hold on;
plot([mu_vec(1),mu_vec(1)],[0,0.6])
xlabel('$p_1$','interpreter','latex')
subplot(1,2,2)
histogram(sample(10000:30000,1),'Normalization','probability')
hold on;
plot([-8,-8],[0,0.1])
xlabel('$log(p_1)$','interpreter','latex')

figure
subplot(1,2,1)
histogram(exp(sample(10000:30000,2)),'Normalization','probability')
hold on;
plot([mu_vec(2),mu_vec(2)],[0,0.6])
xlabel('$p_2$','interpreter','latex')
subplot(1,2,2)
histogram(sample(10000:30000,2),'Normalization','probability')
hold on;
plot([log(mu_vec(2)),log(mu_vec(2))],[0,0.1])
xlabel('$log(p_2)$','interpreter','latex')

figure
subplot(1,2,1)
histogram(exp(sample(10000:50000,3)),'Normalization','probability')
hold on;
plot([t_d,t_d],[0,0.6])
xlabel('$\tau$','interpreter','latex')
subplot(1,2,2)
histogram(sample(10000:50000,3),'Normalization','probability')
hold on;
plot([log(t_d),log(t_d)],[0,0.1])
xlabel('$log(\tau)$','interpreter','latex')

figure
scatterhist(sample(10000:50000,1),sample(10000:50000,2),'kernel','overlay',...
    'Marker','.','MarkerSize',10);
xlabel('$log(p_1)$','interpreter','latex')
ylabel('$log(p_2)$','interpreter','latex')
title('(a) joint posterior of $log(p_1)$ and $log(p_2)$','interpreter','latex')
hold on;

h1 = density2_plot([sample(10000:50000,1),sample(10000:50000,2)],20);
colormap jet;
colorbar
s = scatter(sample(45000:50000,1),sample(45000:50000,2),6,'g','filled');
alpha(s,0.1)
%scatter(log(mu_vec(1)),log(mu_vec(2)),30,'+','r')

figure
scatterhist(sample(10000:50000,1),exp(sample(10000:50000,3)),'kernel','overlay',...
    'Marker','.','MarkerSize',10);
xlabel('$log(p_1)$','interpreter','latex')
ylabel('$\tau$','interpreter','latex')
title('(b) joint posterior of $log(p_1)$ and $\tau$','interpreter','latex')

hold on;
h1 = density2_plot([sample(10000:50000,1),exp(sample(10000:50000,3))],20);
s = scatter(sample(45000:50000,1),exp(sample(45000:50000,3)),6,'g','filled');
alpha(s,0.2)
%scatter(log(mu_vec(1)),t_d,30,'+','r')
colormap jet
colorbar

figure
scatterhist(sample(10000:50000,2),exp(sample(10000:50000,3)),'kernel','overlay',...
    'Marker','.','MarkerSize',10);
xlabel('$log(p_2)$','interpreter','latex')
ylabel('$\tau$','interpreter','latex')
title('(c) joint posterior of $log(p_2)$ and $\tau$','interpreter','latex')
colormap jet;
colorbar
hold on;
h1 = density2_plot([sample(10000:50000,2),exp(sample(10000:50000,3))],20);
sc = scatter(sample(45000:50000,2),exp(sample(45000:50000,3)),6,'g','filled');
alpha(sc,0.2)
%scatter(log(mu_vec(2)),t_d,30,'+','r')