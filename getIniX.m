%% This is the function to obtain design matrix and the corresponding quantity of interest
% Output: logmu: the design matrix
%         Y: the corresponding quantity of interest, here is the square
%         root of observed mutation rate
% Input: a: the parameter for the birth-death process model, for the
%           lifetime of the first generation, which follows a exponential
%           distribution
%        logmu_range: the range to obtain the design matrix
%        chkt: parameter of the birth-death process model, the time
%        checking point to observer the number of cells and the number of
%        mutations
%        num_training: the number of replications at each design points

function [logmu,Y] = getIniX(a,logmu_range,chkt,num_training, num_rep,time_update, L)
    Y = NaN(num_training,1,num_rep);
    logmu = linspace(logmu_range(1),logmu_range(2), num_training)';
    parfor i = 1:num_training
        mu_temp = exp(logmu(i));
        if mod(i, time_update/20) == 0 % screen print to show progress.
           fprintf(['i=', int2str(i),'\n']);
        end
        for j = 1:num_rep
            temp1 = NaN(L,1);
            tmpe2 = NaN(L,2);
            for k = 1:L
              [temp1(k), temp2(k)] = countsizeBDtree(a, mu_temp, chkt);
            end
            Y(i,1,j) = sqrt(sum(temp2)/sum(temp1));
        end
    end
end
function [Nt, Xt] = countsizeBDtree(a, mu, chkt)

Nt = 0;
Xt = 0;
n = 1; % the number of particles in the first generation
lifetime = exprnd(1 / a, [1, n]); % The lifetime of the first generation, following exponential distribution
dtvec = lifetime; % division time or life time
mvec = 0; % is mutant?
f_continue = (dtvec < chkt); % flag of particles that will continue to divide
n_continue = sum(f_continue);
% disp(dtvec);
% disp(mvec);
% disp('--------------------------------------------');
Nt = Nt + sum(~f_continue);
Xt = Xt + sum((~f_continue) & (mvec == 1));
while n_continue > 0
	dtvec_last = dtvec(f_continue);
	mvec_last = mvec(f_continue);
	tmp = [dtvec_last; dtvec_last]; % 2 offsprings
	btvec = tmp(:)';
	lifetime = exprnd(1 / a, [1, 2 * n_continue]);
	dtvec = btvec + lifetime;
	tmp = [mvec_last; mvec_last]; % 2 offsprings
	mvec = binornd(1, (1 - mu) .* tmp(:)' + mu); % mutant always produces mutant offsprings
	% disp(btvec);
	% disp(lifetime);
	% disp(dtvec);
	% disp(mvec);
	% disp('--------------------------------------------');
	f_continue = (dtvec < chkt);
	n_continue = sum(f_continue);
	Nt = Nt + sum(~f_continue);
	Xt = Xt + sum((~f_continue) & (mvec == 1));
end
end
