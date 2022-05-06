function H = simulation1_server(p, Jcase)
  % including folders containing functions to be used
  % The main folder contains the functions to be directly used
  cd('~/ABC/code_rev')
  nSimu = 100;
  POOL = parpool('local',10);

  CILength = NaN(5,3);
  for p = 4:8
    for j = 1:3
      fprintf(['p=', int2str(p),'j=', int2str(j),'\n']);
      l = NaN(nSimu,1);
      u = NaN(nSimu,1);
      for seed = 1:nSimu
         load(strcat('/data/ZChenLab/Lu/ABC/data/simu1/p',num2str(p), ...
          '/Y',num2str(j),'_',num2str(seed),'.mat'))
         [l_tmp, u_tmp] = BootCI_fluc_exp1(Nt, Xt, 0.05, 1000);
         l(seed) = l_tmp;
         u(seed) = u_tmp;
      end
      CILength(p-3,j) = mean(u-l);
    end
  end
  writematrix(CILength,'/data/ZChenLab/Lu/ABC/result/simu1/CILength.csv') 
  
  delete(POOL);
H = 1;
