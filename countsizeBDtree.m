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