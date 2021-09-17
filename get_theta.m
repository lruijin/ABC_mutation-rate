% This is the function to propose new theta 
% Input: 
%   model: 1 constant mutation model; 2. stage-wise mutation rates
%   trans_step: step width of proposal distribution
%   theta_old: theta from last step
%   param_range: the range for parameters
% Output:
%   theta: the new proposed theta
%   param_range_new: the new parameter range

function [theta,param_range_new] = get_theta(model,trans_step,theta_old,param_range)
if model == 1
    pd = makedist('Normal','mu',theta_old,'sigma',trans_step);
    t = truncate(pd,param_range(1), param_range(2));
    theta  = random(t,1,1);
    param_range_new = [min(theta,theta_old)-1e-8,max(theta,theta_old)+1e-8];
elseif model == 2
    p = length(theta_old);
    theta = NaN(1,p);
    param_range_new = NaN(p,2);
    for i = 1:p
        mu_temp = theta_old(i);
        sig_temp = trans_step(i);
        upper_temp = param_range(i,2);
        lower_temp = param_range(i,1);
        if sig_temp == 0
            theta(i) = mu_temp;
            param_range_new(i,:) = param_range(i,:);
        else
            pd = makedist('Normal','mu',mu_temp,'sigma',sig_temp);
            t = truncate(pd,lower_temp, upper_temp);
            theta(i) = random(t,1,1);
            param_range_new(i,:) = [min(theta(i),mu_temp)-1e-8,max(theta(i),mu_temp)+1e-8];
        end
    end
else 
    error('model can only be 1 or 2. 1: constant mutation rate; 2: mutation rate is step funciton');
end
    
end