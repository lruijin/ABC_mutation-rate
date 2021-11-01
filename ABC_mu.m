%% This function is for generating MCMC series of CONSTANT mutation rate based on BD model
% Input: theta_mu: prior for theta, log scale
%        theta_sigma: prior variance for theta, 0 if parameter is not
%                     estimated in this case.
%        param_range; a structure storing range of each parameter
%        obs_X: observed data
%        model_spec: model specifications, including: 
%                     num_training_theta: vector of 2, the number of training data points 
%                                         for initial training and additional trainings.
%                     eps: accuracy threshold, 
%                     psi: error threshold
%                     N: sample size
%                     chkt: time checking point
%                     a: parameter for birth-death process
%                     num_rep: number of replicates at each training point
%                     time_update: for showing progress
%                     init_param: initial hyperparameters for GP surrogates
%                     theta_old: initial values for the parameters to be estimated
%                     bounds: bounds for the parameters to be estimated
%                     trans_step: stepwidth for proposal distribution
function [Sample] = ABC_mu(theta_mu, theta_sigma, param_range,obs_X, model_spec)

S0 = model_spec.num_training_theta(1);
delta = model_spec.num_training_theta(2);
p = 1;   %length(theta_mu)

ksi = model_spec.ksi;
eps = model_spec.eps ;
M = model_spec.M;
N = model_spec.N;
a = model_spec.a;
chkt = model_spec.chkt;

num_rep = model_spec.num_rep;
time_update = model_spec.time_update;
L = model_spec.L;

%kern_param = [0.01,0.4];
init_param = model_spec.init_param;
bounds = model_spec.bounds;


Y_m = mean(obs_X);
Y_v = var(obs_X);

[theta_list, X] = getIniX(a,param_range,chkt,S0, num_rep,time_update,L);

J = size(Y_m,2); % number of the features
[param, model] = mleHomGP(theta_list,X,init_param,bounds);

theta_old = model_spec.theta_old;
trans_step = model_spec.trans_step;
Sample = NaN(N,p);
    
%h = waitbar(0,'Please wait...');
for ss =  1:N
    if mod(ss, time_update) == 0 % screen print to show progress.
       fprintf(['ss=', int2str(ss),'\n']);
    end
    [theta_new,param_range_new] = get_theta(1,trans_step,theta_old,param_range);
    theta = [theta_new;theta_old];
    
    for iter = 1:5000
        alpha = NaN(M,J);  
        [Mu_cond, Var_cond,nugs] = GPHomPrediction(theta,model,param);
        for j = 1:J    
            Yj = Y_m(j);
%             Muj = Mu_cond(:,j);
%             Sigmaj = diag(Var_cond(:,j));
%             sample = randmvn(Muj,Sigmaj,M,eps)';
%             alpha(:,j) = -0.5 * (Yj - sample(:,1)).^2./(nugs(1)/num_rep+eps^2) + 0.5 * (Yj - sample(:,2)).^2./(nugs(2)/num_rep+eps^2);
            sample1 = randn(M,1)*sqrt(Var_cond(1))+repmat(Mu_cond(1),M,1);
            sample2 = randn(M,1)*sqrt(Var_cond(2))+repmat(Mu_cond(2),M,1);
            alpha(:,j) = -0.5 * (Yj - sample1).^2./(nugs(1)/num_rep+eps^2) + 0.5 * (Yj - sample2).^2./(nugs(2)/num_rep+eps^2);
        end
        log_prior_dens = tnorm([theta_new,theta_old], theta_mu, theta_sigma, param_range);
        log_prior_ratio = log_prior_dens(1) - log_prior_dens(2);
        %log_prior_ratio = 0;
        log_like_ratio = sum(alpha,2);
        log_trans_dens_new = tnorm(theta_new, theta_old, trans_step, param_range);
        log_trans_dens_old = tnorm(theta_old, theta_new, trans_step, param_range);
        log_trans_ratio = log_trans_dens_old - log_trans_dens_new; % all cancell out except the jacobian.
        
        Alpha =min([ones(M,1) exp(log_prior_ratio + log_like_ratio + log_trans_ratio)],[],2);
        [Eps,Tau] = unconditionError(Alpha,0.001);
        if Eps <= ksi
            break
        else
            [n_training,~,num_rep] = size(X);
            [Theta_new, X_delta] = getIniX(a,param_range_new,chkt,delta, num_rep,time_update,L);
            [n_training_delta,~,num_rep] = size(X_delta); 
            feat_delta = NaN(n_training_delta + n_training,J,num_rep);
            for k = 1:num_rep
                feat_delta((n_training + 1):(n_training_delta + n_training),:,k) = X_delta(:,:,k);
                feat_delta(1:n_training,:,k) = X(:,:,k);
            end
            X = feat_delta;
            theta_list_new = [theta_list;Theta_new];
            theta_list = theta_list_new;
           
            [param, model] = mleHomGP(theta_list,X,init_param,bounds);
            param_range_new = [param_range_new(1) - 0.001*iter, param_range_new(2)+0.001*iter];
        end
    end
    if rand(1) < Tau
        Sample(ss,:) = theta_new;
    else
        Sample(ss,:) = theta_old;
    end
     theta_old = Sample(ss,:);
     %h = waitbar(ss/N,h,['remaining sample =',num2str(N-ss),'samples']);
end
close(gcf)



