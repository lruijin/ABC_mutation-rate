function [theta_list,Y] = getIniX2(a,param_range,chkt,num_training, num_rep,time_update)
    [theta_list] = getThetaList(num_training,param_range);
    n_theta = size(theta_list,1);
    Y = NaN(n_theta,1,num_rep);
%     if use_parfor == 1
%         parfor i = 1:n_theta
%             theta_temp = theta_list(i,:);
%             mu_vec_temp = [exp(theta_temp(1)) exp(theta_temp(2))];
%             t_temp = exp(theta_temp(3));
%             if mod(i, time_update/20) == 0 % screen print to show progress.
%                fprintf(['i=', int2str(i),'\n']);
%             end
%             for j = 1:num_rep
%                 [temp1, temp2] = countsizeBDtree2(a, mu_vec_temp,t_temp, chkt);
%                 Y(i,1,j) = sqrt(sqrt(temp2/temp1));
%             end
%         end
%     else
        parfor i = 1:n_theta
            theta_temp = theta_list(i,:);
            mu_vec_temp = [exp(theta_temp(1)) exp(theta_temp(2))];
            t_temp = exp(theta_temp(3));
            if mod(i, time_update/20) == 0 % screen print to show progress.
               fprintf(['i=', int2str(i),'\n']);
            end
            for j = 1:num_rep
                [temp1, temp2] = countsizeBDtree2(a, mu_vec_temp,t_temp, chkt);
                Y(i,1,j) = sqrt(sqrt(temp2/temp1));
            end
        end
%     end
end
function[theta_list] = getThetaList(num_training,param_range)
    lhs = lhsdesign(num_training,size(param_range,1));
    theta_list = repmat(param_range(:,1)',num_training,1)+ lhs.*repmat((param_range(:,2) - param_range(:,1))',num_training,1);
end
% function [theta_list] = getThetaList(num_training,param_range)
%     if length(find(num_training>0))==1 % this is the 1-d case
%         n_theta = num_training(num_training > 0); % total number of theta to produce. 
%         theta_list1 = linspace(param_range(1,1), param_range(1,2), n_theta)';
%         theta_list2 = linspace(param_range(2,1), param_range(2,2), n_theta)';
%         theta_list3 = linspace(param_range(3,1), param_range(3,2), n_theta)';
%         theta_list = [theta_list1, theta_list2, theta_list3];
%     else       
%         if  num_training(1)>0
%             theta_list1 = linspace(param_range(1,1), param_range(1,2),num_training(1));  
%         else
%             if param_range(1,1)~=param_range(1,2)
%                 error('param_range does not match with num_training_theta!');
%             end
%             theta_list1 = param_range(1,1); % if 0, we will fix theta at one of its bound.  
%         end
%     
%         if  num_training(2)>0
%             theta_list2 = linspace(param_range(2,1), param_range(2,2),num_training(2));  
%         else
%             if param_range(2,1)~=param_range(2,2)
%                 error('param_range does not match with num_training_theta!');
%             end
%             theta_list2 = param_range(2,1); % if 0, we will fix theta at one of its bound.  
%         end
%     
%         if  num_training(3)>0
%             theta_list3 = linspace(param_range(3,1), param_range(3,2),num_training(3));  
%         else
%             if param_range(3,1)~=param_range(3,2)
%                error('param_range does not match with num_training_theta!');
%             end
%             theta_list3 = param_range(3,1); % if 0, we will fix theta at one of its bound.  
%         end
%     
%     n_theta = prod(num_training(num_training>0));
%     theta_list = NaN(n_theta, 3);
%     tag = 0;
%     for i = 1:length(theta_list1)
%         for j= 1: length(theta_list2)
%             for l = 1:length(theta_list3)
%                 theta_list(tag+1,:) = [theta_list1(i),theta_list2(j),theta_list3(l)];
%                 tag = tag+1;
%             end
%         end
%     end
%     
%     end
% end
