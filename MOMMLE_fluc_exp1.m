function [phat1, phat2] = MOMMLE_fluc_exp1(Z_vec, X_vec)
% Estimate mutation rate for fluctuation experiment with constant mutation rate
% Z_vec: vector of total # of viable cells at tp for J cultures
% X_vec: vector of # of mutants at tp for J cultures
% phat1: MOM/MLE estimator
% phat2: MOM/MLE estimator

Y_vec = Z_vec - X_vec;
Y_bar = mean(Y_vec);
Z_bar = mean(Z_vec);
phat1 = (1 - log(Y_bar) / log(Z_bar)) / 2;
st = max(1e-10, phat1);
fun = @(ph) (1 - 2 * ph) * Y_bar - (1 - ph) * Z_bar ^ (1 - 2 * ph) + ph;
phat2 = fzero(fun, st, optimset('TolX', 1e-10));
end