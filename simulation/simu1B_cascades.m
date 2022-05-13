addpath('/home/xwwu/MyCode');
p_vec = [1e-4, 5e-4, 1e-3, 5e-3, 1e-2];
np = length(p_vec);
J_vec = [10, 20, 30];
nJ = length(J_vec);
nsimu = 100;
nMCMC = 1000;
s = 0.15;
range = [-7.5, -1.5];
ns = 10;
a = 1;
t0 = 11;
nsample = 10;
sigma0 = 0.02;
kparams0 = [1, 1];
eps = 0.005;
% theta_ini = -5; % set on purpose, can start with MOM for better efficiency
burnin = nMCMC * 9 / 10;

clock1 = tic;
for i = 1 : np
    p = p_vec(i);
    for j = 1 : nJ
        J = J_vec(j);
        Z_mat = NaN(nsimu, J);
        X_mat = NaN(nsimu, J);
        phat_MOM_vec = NaN(1, nsimu);
        phat_MLE_vec = NaN(1, nsimu);
        sample_theta_GPS_mat = NaN(nsimu, nMCMC);
        accp_rate_GPS_vec = NaN(1, nsimu);
        total_time_GPS_vec = NaN(1, nsimu);
        train_time_GPS_vec = NaN(1, nsimu);
        phat_GPS_vec = NaN(1, nsimu);
        seed = i * 10 + j;
        rng(seed);
        for k = 1 : nsimu
            [Z_vec, X_vec] = fluc_exp1(a, p, t0, J);
            Z_mat(k, :) = Z_vec;
            X_mat(k, :) = X_vec;

            [phat_MOM, phat_MLE] = MOMMLE_fluc_exp1(Z_vec, X_vec);
            phat_MOM_vec(k) = phat_MOM;
            phat_MLE_vec(k) = phat_MLE;
            theta_ini = log10(phat_MOM);
            clock2 = tic;
            [sample_theta_GPS, accp_rate_GPS, train_time_GPS] = ABC_fluc_exp1(nMCMC, Z_vec, X_vec, theta_ini, s, range, ns, a, t0, J, nsample, sigma0, kparams0, eps, true);
            total_time_GPS = toc(clock2);
            phat_GPS = 10 ^ mean(sample_theta_GPS((burnin + 1) : end));
            sample_theta_GPS_mat(k, :) = sample_theta_GPS;
            accp_rate_GPS_vec(k) = accp_rate_GPS;
            train_time_GPS_vec(k) = train_time_GPS;
            total_time_GPS_vec(k) = total_time_GPS;
            phat_GPS_vec(k) = phat_GPS;
        end
        save(strcat('/home/xwwu/MyResult/simu1B_', int2str(seed), '.mat'));
    end
end
toc(clock1)
