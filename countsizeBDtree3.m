%% Function of simulating population size and number of mutants
%% This function allows different growth rates and different mutation rates
%% Input: a: growth rate of normal cells
%%        a_dt: delta of a, growth rate of mutants: a + delta
%%        mu_par: parameters concerning mutation rates, could be only one or
%%                three paramters in the stagewise function
%%        chkt: time checking point for counting mutants.
%% Output: Nt: population size at time t
%%         Xt: number of mutants.

function [Nt, Xt] = countsizeBDtree3(a, a_dt, mu_par, chkt)
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
  if length(mu_par) == 1
    mu_vec = [mu_par mu_par];
    t_d = 0;
  elseif length(mu_par) == 3
    mu_vec = [mu_par(1) mu_par(2)];
    t_d = mu_par(3);
  end
 
  while n_continue > 0
    dtvec_last = dtvec(f_continue);
    mvec_last = mvec(f_continue);
    tmp = [dtvec_last; dtvec_last]; % 2 offsprings
    btvec = tmp(:)';
    tmp = [mvec_last; mvec_last];
    lifetime = exprnd(1 ./ (a * (a_dt .^ tmp(:)')));
	  dtvec = btvec + lifetime;
    mu_f = (dtvec > t_d) + 1;% flag showing which mu should be followed
    mu = mu_vec(mu_f);
	  mvec = binornd(1, (1 - mu) .* tmp(:)' + mu); % mutant always produces mutant offsprings
    % disp(btvec);
    % disp(lifetime);
    % disp(dtvec);
    % disp(mvec);
    % :wqdisp('--------------------------------------------');
    f_continue = (dtvec < chkt);
    n_continue = sum(f_continue);
    Nt = Nt + sum(~f_continue);
    Xt = Xt + sum((~f_continue) & (mvec == 1));
  end
  