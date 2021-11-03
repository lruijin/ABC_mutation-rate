function [Z, X] = mut_bMBP(a, p, t0)
% Generate (z, x) data for bMBP model with constant mutation
% a: rate parameter of exponential life time
% p: mutation probability of each single particle
% t0: time of plating
% Z: total # of viable cells at t0
% X: # of mutants at t0

Z = 0;
X = 0;
dtvec = exprnd(1 / a, [1, 1]); % unit initial size, death time
mvec = 0; % is mutant?
f_continue = (dtvec < t0); % flag of particles that will continue to divide
n_continue = sum(f_continue);
Z = Z + sum(~f_continue);
X = X + sum((~f_continue) & (mvec == 1));
while n_continue > 0
	dtvec_last = dtvec(f_continue);
	mvec_last = mvec(f_continue);
	dtvec = repelem(dtvec_last, 2) + exprnd(1 / a, [1, 2 * n_continue]); % 2 offsprings
% 	dtvec = reshape([dtvec_last; dtvec_last], 1, []) + exprnd(1 / a, [1, 2 * n_continue]); % for old version Matlab
	mvec = binornd(1, (1 - p) .* repelem(mvec_last, 2) + p); % mutant always produces mutant offsprings
% 	mvec = binornd(1, (1 - p) .* reshape([mvec_last; mvec_last], 1, []) + p); % for old version Matlab
	f_continue = (dtvec < t0);
	n_continue = sum(f_continue);
	Z = Z + sum(~f_continue);
	X = X + sum((~f_continue) & (mvec == 1));
end
end