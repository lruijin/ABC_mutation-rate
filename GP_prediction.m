function [Mu_cond, Var_cond, Kern2] = GP_prediction(kern_param, true_theta, theta_list, init_feat_s, Kern1)

% This function predict the feature value at theta_true based 
% on existing theta_list, and init_feat. 
% Input:
% kern_param: a 1 by 2 vector. Hyperparameters for the kernel function. See
%            se_kern_fast for definition. 
% theta_true: a n by 3 vector (leaf density, mean diameter, orientation).
%             n is the number of true theta considered.
% theta_list: the initial grid of theta on which to sample initial echoes.
% ini_feature: the intial feature extracted. 
% rescale_theta: whether or not to rescale theta to [0,1].
% feat_idx: if an integer, only consider prediction for that index of the
%           feature vector. if empty [], we will do all features one by
%           one (assuming features are independent of each other).
% Kern1: the Koldold value that have been precalculated. 

% Output:
% Mu_cond: a n by length(feat_idx) vector,
%          n is the number of true theta considered
% Var_cond: a n by n by lenghth(feat_idx) array,
%           n is the number of true tehta considered
% Kern: the structure of storing the kernel values: K_new_new, K_new_old,
%       K_old_old
% Reference:
% Chapter 2: Rasmussen and Williams Gaussian Processes for Machine Learning

% step 2: calculate kernel, K(TH, TH), K(TH, th), K(th, th)
% here I assume th is a scalar, on which I want to predict f*. 
% OK, I think I will treat th as a vector, that will save me a lot of
% computation time. 
% sig2_n= min(1e-4, var(init_feat_s)*0.01); 
% sig2_n = min(1e-6,var(y)); % this parameter is extremely important, it influences the confidence interval of predicted feature values.

% true_feat=wave.D(:,wave.comp_list);
% true_feat1=true_feat(:,feat_idx);
% sig2_n = min(1e-4,var(true_feat1)*0.0005);

n = size(true_theta,1);

K_new_old = se_kern_fast(kern_param, true_theta, theta_list);
K_new_new = se_kern_fast(kern_param, true_theta);

% Mu_cond = NaN(1, length(feat_idx));
% Var_cond = NaN(1, length(feat_idx));
J = size(init_feat_s,2);
Mu_cond = NaN(n, J);
Var_cond = NaN(n, n, J);
K_iv = Kern1.K_iv;
% if use_parfor ==1
%     parfor i = 1:J
%         K_ivUoo = K_iv(:,:,i);
%         Mu_cond(:,i) = K_new_old * (K_ivUoo*(K_ivUoo'*init_feat_s(:,i))); 
%         Var_cond(:,:,i) =  K_new_new - (K_new_old*K_ivUoo)*(K_ivUoo'*K_new_old'); 
%     end
% else
    for i = 1:J
        K_ivUoo = K_iv(:,:,i);
        Mu_cond(:,i) = K_new_old * (K_ivUoo*(K_ivUoo'*init_feat_s(:,i))); 
        Var_cond(:,:,i) =  K_new_new - (K_new_old*K_ivUoo)*(K_ivUoo'*K_new_old'); 
    end
%end
Kern2.new_old = K_new_old;
Kern2.old_old = Kern1.old_old;
Kern2.new_new = K_new_new;
 
 