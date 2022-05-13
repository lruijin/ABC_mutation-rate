function [l, u] = BootCI_fluc_exp1(Z_vec, X_vec, alpha, nboot)

J = length(Z_vec);
s = RandStream('mlfg6331_64');
phat_vec = NaN(1, nboot);
parfor i = 1 : nboot
    ix = randsample(s, J, J, true);
    Z_boot = Z_vec(ix);
    X_boot = X_vec(ix);
    [~, phat] = MOMMLE_fluc_exp1(Z_boot, X_boot);
    phat_vec(i) = phat;
end
l = quantile(phat_vec, alpha / 2);
u = quantile(phat_vec, 1 - alpha / 2);
end