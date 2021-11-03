% For Figure 1, choose summary statistic for parallel cultured data with constant mutation
addpath('C:/Users/Xiaowei/Documents/Work/MutationProject/MyCode');
a = 1;
p_vec = 1e-5 : 8e-5 : 1e-2;
np = length(p_vec);
t0 = 11;
J = 30;
S1_vec = NaN(1, np); % [1 - log(bar{Y})/log(bar{Z})] / 2
S2_vec = NaN(1, np); % log(Y/Z)
S3_vec = NaN(1, np); % sqrt(X/Z)
rng(1);
tic;
for i = 1 : np
    p = p_vec(i);
    [Z_vec, X_vec] = fluc_exp1(a, p, t0, J);
    S1_vec(i) = (1 - log(mean(Z_vec - X_vec)) / log(mean(Z_vec))) / 2;
    S2_vec(i) = mean(log(1 - X_vec ./ Z_vec));
    S3_vec(i) = mean(sqrt(X_vec ./ Z_vec));
end
toc;
figure;
subplot(2, 2, 1);
plot(p_vec, S1_vec, 'k-');
title('(A) $\frac{1}{2}[1 - \log(\overline{Y})/\log(\overline{Z})]$ VS. $p$', 'interpreter', 'latex');
hold on;
plot([min(p_vec), max(p_vec)], [min(p_vec), max(p_vec)], 'k--');
subplot(2, 2, 2);
plot(p_vec, S2_vec, 'k-');
title('(B) $\overline{\log(Y/Z)}$ VS. $p$', 'interpreter', 'latex');
subplot(2, 2, 3);
plot(p_vec, S3_vec, 'k-');
title('(C) $\overline{\sqrt{X/Z}}$ VS. $p$', 'interpreter', 'latex');
subplot(2, 2, 4);
plot(log10(p_vec), S3_vec, 'k-');
title('(D) $\overline{\sqrt{X/Z}}$ VS. $\log_{10}(p)$', 'interpreter', 'latex');
set(gcf, 'PaperUnits', 'centimeters', 'PaperPosition', [0 0 22 12]);
saveas(gcf, 'C:/Users/Xiaowei/Documents/Work/MutationProject/MyResult/Fig1', 'epsc');
