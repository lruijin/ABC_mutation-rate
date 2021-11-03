% run on server, J = 10 takes 3e4 seconds = 9.38 hrs
% nohup matlab -nodisplay -nosplash -r "seed = 11" <./simu1A_cascades.m> ./temp11.txt &
% nohup matlab -nodisplay -nosplash -r "seed = 12" <./simu1A_cascades.m> ./temp12.txt &
% nohup matlab -nodisplay -nosplash -r "seed = 13" <./simu1A_cascades.m> ./temp13.txt &
% nohup matlab -nodisplay -nosplash -r "seed = 21" <./simu1A_cascades.m> ./temp21.txt &
% nohup matlab -nodisplay -nosplash -r "seed = 22" <./simu1A_cascades.m> ./temp22.txt &
% nohup matlab -nodisplay -nosplash -r "seed = 23" <./simu1A_cascades.m> ./temp23.txt &
% nohup matlab -nodisplay -nosplash -r "seed = 31" <./simu1A_cascades.m> ./temp31.txt &
% nohup matlab -nodisplay -nosplash -r "seed = 32" <./simu1A_cascades.m> ./temp32.txt &
% nohup matlab -nodisplay -nosplash -r "seed = 33" <./simu1A_cascades.m> ./temp33.txt &
% nohup matlab -nodisplay -nosplash -r "seed = 41" <./simu1A_cascades.m> ./temp41.txt &
% nohup matlab -nodisplay -nosplash -r "seed = 42" <./simu1A_cascades.m> ./temp42.txt &
% nohup matlab -nodisplay -nosplash -r "seed = 43" <./simu1A_cascades.m> ./temp43.txt &
% nohup matlab -nodisplay -nosplash -r "seed = 51" <./simu1A_cascades.m> ./temp51.txt &
% nohup matlab -nodisplay -nosplash -r "seed = 52" <./simu1A_cascades.m> ./temp52.txt &
% nohup matlab -nodisplay -nosplash -r "seed = 53" <./simu1A_cascades.m> ./temp53.txt &

addpath('/home/xwwu/MyCode');
parpool(8);
load(strcat('/home/xwwu/MyResult/simu1B_', int2str(seed), '.mat'));
sample_theta_ABC_mat = NaN(nsimu, nMCMC);
accp_rate_ABC_vec = NaN(1, nsimu);
total_time_ABC_vec = NaN(1, nsimu);
phat_ABC_vec = NaN(1, nsimu);

clock1 = tic;
parfor i = 1 : nsimu
    Z_vec = Z_mat(i, :);
    X_vec = X_mat(i, :);
    theta_ini = log10(phat_MOM_vec(i));

    clock2 = tic;
    [sample_theta_ABC, accp_rate_ABC, ~] = ABC_fluc_exp1(nMCMC, Z_vec, X_vec, theta_ini, s, range, ns, a, t0, J, nsample, sigma0, kparams0, eps, false);
    total_time_ABC = toc(clock2);
    phat_ABC = 10 ^ mean(sample_theta_ABC((burnin + 1) : end));
    sample_theta_ABC_mat(i, :) = sample_theta_ABC;
    accp_rate_ABC_vec(i) = accp_rate_ABC;
    total_time_ABC_vec(i) = total_time_ABC;
    phat_ABC_vec(i) = phat_ABC;
end
toc(clock1);

save(strcat('/home/xwwu/MyResult/simu1A_', int2str(seed), '.mat'));
delete(gcp('nocreate'));
