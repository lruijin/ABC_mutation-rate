% This function is for getting training theta list as the input for bMBP model with piece-wise constant mutation rates.
% Input: a: rate parameter of the bMBP model.
%        param_range; a structure storing range of each parameter
%        chkt: time checking point
%        num_training: vector of 2, the number of training data points for initial training and additional trainings. 
%        num_rep: number of replicates at each training point             
%        time_update: for showing progress 
%        L: number of bacteria in each culture
% Output: theta_list: training points
%         Y: Output from model

function [theta_list,Y] = getIniX2(a, dt_range,p_range,chkt,num_training, num_rep,time_update, L, constraint)
    if length(dt_range) == 1
        a_dt = dt_range;
        dt_par = 0;
        param_range = p_range;
    else 
        dt_par = 1;
        param_range = [dt_range;p_range];
    end
    
    p_par = size(p_range,1);
    
    [theta_list] = getThetaList(num_training, param_range, dt_par, constraint);
    n_theta = size(theta_list,1);
    Y = NaN(n_theta,1,num_rep);
    parfor i = 1:n_theta
        theta_temp = theta_list(i,:);
        if dt_par == 0
            if p_par == 1
                mu_par_temp = 10^theta_temp(1);
            else
                mu_par_temp = [10^(theta_temp(1)) 10^(theta_temp(2)) theta_temp(3)];
            end
        else 
            a_dt = 10^theta_temp(1);
            if p_par == 1;
                mu_par_temp = 10^(theta_temp(2));
            else
                mu_vec_temp = [10^(theta_temp(2)) 10^(theta_temp(3))];
                t_temp = theta_temp(4);
                mu_par_temp =[mu_vec_temp t_temp];
            end
        end
                
        if mod(i, time_update/50) == 0 % screen print to show progress.
            fprintf(['i=', int2str(i),'\n']);
        end
        for j = 1:num_rep
            tmp1 = NaN(L,1);
            tmp2 = NaN(L,1);
            for k = 1:L
                [tmp1(k), tmp2(k)] = countsizeBDtree3(a, a_dt, mu_par_temp, chkt);
            end
            Y(i,1,j) = sqrt(sqrt(sum(tmp2)/sum(tmp1)));
        end
    end
        
%     end
end

function[theta_list] = getThetaList(num_training, param_range, dt_par, constraint)
    lhs = lhsdesign(num_training,size(param_range,1));
    theta_list = repmat(param_range(:,1)',num_training,1) + ...
            lhs.*repmat((param_range(:,2) - param_range(:,1))',num_training,1);
    if constraint == 1
        if dt_par == 0
            idx = theta_list(:,2) > theta_list(:,1);
        else
            idx = theta_list(:,3) > theta_list(:,2);
        end
    end
    
end
