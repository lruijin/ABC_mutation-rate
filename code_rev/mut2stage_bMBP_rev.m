function [Z, X] = mut2stage_bMBP_rev(Z0, a, delta, p1, p2, tj, tp)
% Generate (z, x) data for bMBP model with constant mutation
% Z0: # of non-mutants at t = 0
% a: rate parameter of exponential life time for non-mutants
% delta: growth parameter for mutants relative to non-mutants
% p1: mutation probability in stage 1
% p2: mutation probability in stage 2
% tj: jumping time of the piece-wise constant mutation probability
% tp: time of plating
% Z: total # of viable cells at tp
% X: # of mutants at tp

[Z1, X1] = mut_bMBP_rev(Z0, a, delta, p1, tj);
[Z, X2] = mut_bMBP_rev(Z1 - X1, a, delta, p2, tp - tj);
X = sum(geornd(exp(-(a * delta) .* (tp - tj)), 1, X1)) + X1 + X2;
end