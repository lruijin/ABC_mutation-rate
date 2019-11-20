cd('C:\Users\lruijin\Documents\MATLAB\ABC_wu\Code')
addpath('C:\Users\lruijin\Documents\MATLAB\ABC_wu\Code\kernel')
theta_test = linspace(1e-4,0.009,100);
logtheta_test = linspace(-10,-2,100);

chkt = 9; a = 1;
Nt1 = NaN(1000,100);
Xt1 = NaN(1000,100);
Nt2 = NaN(1000,100);
Xt2 = NaN(1000,100);
for i = 1:100
    theta = theta_test(i);
    Yt = NaN(1000,1);
    logYt = NaN(1000,1);
    for j = 1:1000
        [temp1, temp2] = countsizeBDtree(a, theta, chkt);%%
        Nt1(j,i) = temp1;
        Xt1(j,i) = temp2;
    end
    fprintf(['i=', int2str(i),'\n']);
end

for i = 1:100
    theta = exp(logtheta_test(i));
    for j = 1:1000
        [temp1, temp2] = countsizeBDtree(a, theta, chkt);%%
        Nt2(j,i) = temp1;
        Xt2(j,i) = temp2;
    end
    fprintf(['i=', int2str(i),'\n']);
end

X1 = mean(Xt1,1);
logX1 = mean(log(Nt1 - Xt1),1);
R1 = mean(Xt1./Nt1,1);
logR1 = mean(log(1-Xt1./Nt1),1);
square_root1 = mean(sqrt(Xt1./Nt1),1);
mom1 = mean(log(Nt1-Xt1)./log(Nt1),1);

X2 = mean(Xt2,1);
logX2 = mean(log(Nt2-Xt2),1);
R2 = mean(Xt2./Nt2,1);
logR2 = mean(log(1-Xt2./Nt2),1);
square_root2 = mean(sqrt(Xt2./Nt2),1);
mom2 = mean(log(Nt2-Xt2)./log(Nt2),1);

figure
subplot(2,3,1)
plot(logtheta_test,mom2);
ylabel('$\frac{log(N_t-Z_t)}{log(N_t)}$','interpreter','latex');
xlabel('log(p)');
title('1(a) MOM vs. log(p)','interpreter','latex')
% subplot(2,2,2)
% plot(logtheta_test,logX2);
% ylabel('log(# of mutants)');xlabel('log(mutation rate)');
% subplot(2,2,3)
% plot(logtheta_test,R2);
% ylabel('percent of mutants'); xlabel('log(mutation rate)');
subplot(2,3,2)
plot(logtheta_test,logR2);
ylabel('$log(\frac{N_t-Z_t}{N_t})$','interpreter','latex');
xlabel('log(p)');
title('1(b) $log(\frac{N_t-Z_t}{N_t})$ vs. log(p)','interpreter','latex')

subplot(2,3,3)
plot(logtheta_test,square_root2);
ylabel('$\sqrt{Z_t/N_t}$','interpreter','latex');
xlabel('log(p)');
title('1(c) $\sqrt{Z_t/N_t}$ vs. log(p)','interpreter','latex')

subplot(2,3,4)
plot(theta_test,mom1);
ylabel('$\frac{log(N_t-Z_t)}{log(N_t)}$','interpreter','latex');
xlabel('p');
title('2(a) MOM vs. log(p)','interpreter','latex')

subplot(2,3,5)
plot(theta_test,logR1);
ylabel('$log(\frac{N_t-Z_t}{N_t})$','interpreter','latex');
xlabel('p');
title('2(b) $log(\frac{N_t-Z_t}{N_t})$ vs. p','interpreter','latex')

subplot(2,3,6)
plot(theta_test,square_root1);
ylabel('$\sqrt{Z_t/N_t}$','interpreter','latex');
xlabel('p');
title('2(c) $\sqrt{Z_t/N_t} vs. p$','interpreter','latex')

%% plot the scatter plot of data
figure;
plot(logtheta_test,square_root2);
hold on;
for i=1:length(logtheta_test)
    if(mod(i,5) ==0 )
     plot(logtheta_test(i),sqrt(Xt2(1:50,i)./Nt2(1:50,i)),'.')
    end
end
hold off;
ylabel('$\sqrt{percent\ of\ mutants}$','interpreter','latex');
xlabel('log(mutation rate)');

%% plot the prediction band and the data
save('..\Result\simulation.mat','Nt1','Nt2','Xt1','Xt2','X1','logX1','R1','logR1','X2','logX2','R2','logR2','theta_test','logtheta_test')

S0 = 50;
param_range = [-20,-12];
chkt = 9; a = 1;num_rep = 10;time_update = 5;

[theta_list, X] = getIniX(a,param_range,chkt,S0, num_rep,time_update);
% Z0 = mean(Z,3);
%     mult = size(Z,3);
%     n = size(theta_list,1);
%     Z = reshape(Z,n,mult);
%     Z = reshape(Z',numel(Z),1);
% save('../Result/GP_data.mat','theta_list','X','Z','Z0');
init_param = [1;1e-4];
bounds.upper = 8;
bounds.lower = 1e-4;

[param, model] = mleHomGP(theta_list,mean(X,3),init_param,bounds);
load('../Result/simulation.mat')
[est_mu, est_var,nugs] = GPHomPrediction(theta_list,model,param);

figure
mycolor = jet;

plot(logtheta_test,mean(sqrt(Xt2./Nt2)),'color',mycolor(10,:));
hold on
plot(logtheta_test,est_mu,'color',mycolor(50,:));
plot(logtheta_test,est_mu+norminv(0.025,0,sqrt(est_var)),'--','color',mycolor(50,:));
plot(logtheta_test,est_mu+norminv(0.975,0,sqrt(est_var)),'--','color',mycolor(50,:));
for i=1:5
     dot = scatter(logtheta_test,mean(sqrt(Xt2((20*(i-1)+1):(20*i),:)./Nt2((20*(i-1)+1):(20*i),:))),'filled','sizeData',10);
     alpha(dot,0.3)
end
ylabel('$\sqrt{percent\ of\ mutants}$','interpreter','latex');
xlabel('log(p)');

lgd = legend('sample mean','prediction mean','95% predictive interval');
lgd.FontSize=8;
set(lgd,'Location','northwest')

hold off

%% --------------------Case 2: for the stage function model----------------
S0 = 150;
param_range = [-9,-1;-9,-1;-0.1,2.2];
a =  1; chkt =9;num_rep = 10;time_update = 50;
[theta_list,Y] = getIniX2(a,param_range,chkt,S0, num_rep,time_update);

hFig = figure(); %set(gcf,'color','w');
winsize = get(hFig,'Position');
winsize1 = [10, 500, winsize(3)*2, winsize(3)*0.8];
set(hFig, 'Position', winsize1);


mycolor = jet;
a1 = subplot(1,3,1);
set(a1, 'position', [0.1, 0.200, 0.25, 0.71])
scatter3(theta_list(:,1),theta_list(:,2),exp(theta_list(:,3)),30*mean(Y,3),mycolor(round(size(mycolor,1)*mean(Y,3)),:),'filled');
xlabel('$log(p_1)$','interpreter','latex');
ylabel('$log(p_2)$','interpreter','latex');
zlabel('$\tau$','interpreter','latex');
colormap(mycolor);
col = colorbar;
set(col,'position',[0.04, 0.200, 0.01, 0.71]);
col.Label.String = '$\sqrt[4]{mutants \ proportion}$';
col.Label.Interpreter = 'latex';

a2 = subplot(1,3,2);
set(a2, 'position', [0.4, 0.200, 0.25, 0.71])
scatter3(theta_list(:,1),theta_list(:,2),exp(theta_list(:,3)),30*mean(Y,3),mycolor(round(size(mycolor,1)*mean(Y,3)),:),'filled');
xlabel('$log(p_1)$','interpreter','latex');
ylabel('$log(p_2)$','interpreter','latex');
zlabel('$\tau$','interpreter','latex');

a3 = subplot(1,3,3);
set(a3, 'position', [0.7, 0.200, 0.25, 0.71])
al3 = scatter3(theta_list(:,1),theta_list(:,2),exp(theta_list(:,3)),30*mean(Y,3),mycolor(round(size(mycolor,1)*mean(Y,3)),:),'filled');
xlabel('$log(p_1)$','interpreter','latex');
ylabel('$log(p_2)$','interpreter','latex');
zlabel('$\tau$','interpreter','latex');

% GP prediction based on the 150 training points, 10 replications each
init_param = [1; 1; 1; 1e-4];
bounds.upper = [8 8 8];
bounds.lower = [1e-4 1e-4 1e-4];
[param, model] = mleHomGP(theta_list,Y,init_param,bounds);

[theta_list0,Y0] = getIniX2(a,param_range,chkt,S0, 500,time_update);
[mu, var,nugs] = GPHomPrediction(theta_list0,model,param);

hFig = figure(); %set(gcf,'color','w');
winsize = get(hFig,'Position');
winsize1 = [10, 500, winsize(3)*2, winsize(3)*0.8];
set(hFig, 'Position', winsize1);

a1 = subplot(1,3,1);
set(a1, 'position', [0.1, 0.200, 0.25, 0.71])
plot3d_errorbars(theta_list0(:,1), theta_list0(:,2), mu, 1.96*sqrt(var));
scatter3(theta_list0(:,1), theta_list0(:,2), mean(Y0,3),15,'red','filled');
xlabel('$log(p_1)$','interpreter','latex');
ylabel('$log(p_2)$','interpreter','latex');
zlabel('$\sqrt[4]{mutants \ proportion}$','interpreter','latex');

a2 = subplot(1,3,2);
set(a2, 'position', [0.4, 0.200, 0.25, 0.71])
plot3d_errorbars(theta_list0(:,1), exp(theta_list0(:,3)), mu, 1.96*sqrt(var));
scatter3(theta_list0(:,1), exp(theta_list0(:,3)), mean(Y0,3),15,'red','filled');
xlabel('$log(p_1)$','interpreter','latex');
ylabel('$\tau$','interpreter','latex');
zlabel('$\sqrt[4]{mutants \ proportion}$','interpreter','latex');

a3 = subplot(1,3,3);
set(a3, 'position', [0.7, 0.200, 0.25, 0.71])
plot3d_errorbars(theta_list0(:,2), exp(theta_list0(:,3)), mu, 1.96*sqrt(var));
scatter3(theta_list0(:,2), exp(theta_list0(:,3)), mean(Y0,3),15,'red','filled');
xlabel('$log(p_2)$','interpreter','latex');
ylabel('$\tau$','interpreter','latex');
zlabel('$\sqrt[4]{mutants \ proportion}$','interpreter','latex');