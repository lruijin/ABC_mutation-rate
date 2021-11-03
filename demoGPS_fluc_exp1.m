% For Figure 2, demonstrate GP regression for multiple samples of parallel cultured data with constant mutation
addpath('C:/Users/Xiaowei/Documents/Work/MutationProject/MyCode');
a = 1;
theta_vec = (-5 : 0.03 : -2)';
p_vec = 10 .^ theta_vec;
t0 = 11;
J = 30;
nsample = 10;
sigma0 = 0.02;
kparams0 = [1, 1];
rng(1);
tic;
[gprMd, theta_rep_vec, S_vec] = trainGPS(a, p_vec, t0, J, nsample, sigma0, kparams0);
toc;
[Spred, Ssd, Sint] = predict(gprMd, theta_vec);
figure();
plot(theta_rep_vec, S_vec, 'r.', 'MarkerSize', 12);
hold on
plot(theta_vec, Spred, 'b-', 'LineWidth', 1.5);
h = fill([theta_vec', fliplr(theta_vec')], [Sint(:, 1)', fliplr(Sint(:, 2)')], [0.8, 0.8, 0.8], 'LineStyle', '--');
set(h, 'facealpha', 0.5);
xlabel('Mutation rate in $\log_{10}$ scale', 'interpreter', 'latex');
ylabel('Summary statistic VS. Predictive mean');
legend({'$\overline{\sqrt{x/z}}$, 10 replicates at each grid point', 'Predictive mean in GP regression', '95\% predictive interval in GP regression'}, 'interpreter', 'latex', 'Location', 'Best');
hold off;
set(gcf, 'PaperUnits', 'centimeters', 'PaperPosition', [0 0 16 12]);
saveas(gcf, 'C:/Users/Xiaowei/Documents/Work/MutationProject/MyResult/Fig2', 'epsc');
