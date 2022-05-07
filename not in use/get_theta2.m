% This is the function for proposing new theta 
% Input: 
  %   model: 1. constant mutation model; 2. piece-wise constant mutation rates
%   trans_step: step width of proposal distribution
%   theta_old: theta from last step
%   param_range: the range for parameters
% Output:
  %   theta: the new proposed theta
%   param_range_new: the new parameter range

function [theta, dt_range_new, p_range_new] = get_theta2(trans_step,theta_old,dt_range,p_range,constraint)
  p_par = size(p_range,1);
  if length(dt_range) == 1
    dt_range_new = dt_range;
    if p_par == 1
      pd = makedist('Normal','mu',theta_old,'sigma',trans_step);
      t = truncate(pd,p_range(1), p_range(2));
      theta  = random(t,1,1);
      p_range_new = [min(theta,theta_old)-1e-8,max(theta,theta_old)+1e-8];
    elseif p_par == 3
      theta = [1 0 0];
      p_range_new = NaN(p_par,2);
      if constraint == 1
        while theta(2) > theta(1)
          for i = 1:2
            mu_temp = theta_old(i);
            sig_temp = trans_step(i);
            upper_temp = p_range(i,2);
            lower_temp = p_range(i,1);
            pd = makedist('Normal','mu',mu_temp,'sigma',sig_temp);
            t = truncate(pd,lower_temp, upper_temp);
            theta(i) = random(t,1,1);
            p_range_new(i,:) = [min(theta(i),mu_temp)-1e-8,max(theta(i),mu_temp)+1e-8];
          end
        end
        mu_temp = theta_old(3);
        sig_temp = trans_step(3);
        upper_temp = p_range(3,2);
        lower_temp = p_range(3,1);
        pd = makedist('Normal','mu',mu_temp,'sigma',sig_temp);
        t = truncate(pd,lower_temp, upper_temp);
        theta(3) = random(t,1,1);
        p_range_new(3,:) = [min(theta(3),mu_temp)-1e-8,max(theta(3),mu_temp)+1e-8];
      elseif constraint == 0
        for i = 1:p_par
          mu_temp = theta_old(i);
          sig_temp = trans_step(i);
          upper_temp = p_range(i,2);
          lower_temp = p_range(i,1);
          pd = makedist('Normal','mu',mu_temp,'sigma',sig_temp);
          t = truncate(pd,lower_temp, upper_temp);
          theta(i) = random(t,1,1);
          p_range_new(i,:) = [min(theta(i),mu_temp)-1e-8,max(theta(i),mu_temp)+1e-8];
        end
      else
        error('constraint can only be 0 or 1. 1: mu1 smaller than mu2; 0: no constraint');
      end
    else
      error("The number of parmaters could be either 1 or 3. 1: single mutation rate; 3:stagewise function for mutation rate");
    end
  elseif length(dt_range) == 2
    pd = makedist('Normal','mu',theta_old(1),'sigma',trans_step(1));
    t = truncate(pd,dt_range(1), dt_range(2));
    theta  = random(t,1,1);
    dt_range_new = [min(theta,theta_old(1))-1e-8,max(theta,theta_old(1))+1e-8];
    if p_par == 1
      pd = makedist('Normal','mu',theta_old(2),'sigma',trans_step(2));
      t = truncate(pd,p_range(1), p_range(2));
      theta_tmp = random(t,1,1);
      p_range_new = [min(theta_tmp,theta_old(2))-1e-8,max(theta_tmp,theta_old(2))+1e-8];
      theta = [theta theta_tmp];
    elseif p_par == 3
      theta_tmp = [1 0 0];
      p_range_new = NaN(p_par,2);
      if constraint == 1
        while theta_tmp(2) > theta_tmp(1)
          for i = 1:2
            mu_temp = theta_old(i+1);
            sig_temp = trans_step(i+1);
            upper_temp = p_range(i,2);
            lower_temp = p_range(i,1);
            pd = makedist('Normal','mu',mu_temp,'sigma',sig_temp);
            t = truncate(pd,lower_temp, upper_temp);
            theta_tmp(i) = random(t,1,1);
            p_range_new(i,:) = [min(theta_tmp(i),mu_temp)-1e-8,max(theta_tmp(i),mu_temp)+1e-8];
          end
        end
        mu_temp = theta_old(4);
        sig_temp = trans_step(4);
        upper_temp = p_range(3,2);
        lower_temp = p_range(3,1);
        pd = makedist('Normal','mu',mu_temp,'sigma',sig_temp);
        t = truncate(pd,lower_temp, upper_temp);
        theta_tmp(3) = random(t,1,1);
        p_range_new(3,:) = [min(theta_tmp(3),mu_temp)-1e-8,max(theta_tmp(3),mu_temp)+1e-8];
        theta = [theta theta_tmp];
      elseif constraint == 0
        for i = 1:p_par
          mu_temp = theta_old(i+1);
          sig_temp = trans_step(i+1);
          upper_temp = p_range(i,2);
          lower_temp = p_range(i,1);
          pd = makedist('Normal','mu',mu_temp,'sigma',sig_temp);
          t = truncate(pd,lower_temp, upper_temp);
          theta_tmp(i) = random(t,1,1);
          p_range_new(i,:) = [min(theta_tmp(i),mu_temp)-1e-8,max(theta_tmp(i),mu_temp)+1e-8];
        end
        theta = [theta theta_tmp];
      else
        error('constraint can only be 0 or 1. 1: mu1 smaller than mu2; 0: no constraint');
      end
    else
      error("The number of mutation rate parmaters could be either 1 or 3. 1: single mutation rate; 3:stagewise function for mutation rate")
    end
  else
    error("The range of delta mean life could either be a number or a vector of 2 numbers.");
  end
end