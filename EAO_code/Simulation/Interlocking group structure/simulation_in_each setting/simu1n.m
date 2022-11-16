  B = 100;

  g1 = 400;
  r1=0.2;
  d1=20;
  inf =  zeros(B,5);
  solve_time = zeros(B,3);
  run_time = zeros(B,3);
  RSS = zeros(B,3);
   parpool('local',100)
  parfor b = 1:B 
              seed = 159866;
          seed = round(0.626*seed - b);
          n = 200;
          tmprun = [b,n, d1, g1, r1,seed];
          disp(tmprun)
          res_cv = simures(n, d1, g1, r1,seed);
          
          inf(b,:) = [n, d1, g1, r1,seed];
          solve_time(b,:) = res_cv(1,:);
          run_time(b,:) = res_cv(2,:);
          RSS(b,:) = res_cv(3,:);
      
  end
  resn200 = [inf, solve_time, run_time, RSS];
  
  csvwrite("newsimudata/resn200_1.csv", [inf, solve_time, run_time, RSS])
  
  fprintf("n one is done")
  
  inf =  zeros(B,5);
  solve_time = zeros(B,3);
  run_time = zeros(B,3);
  RSS = zeros(B,3);
  
  parfor b = 1:B
              seed = 5235;
          seed = round(5.147*seed - b);
          n = 400;
          tmprun = [b,n, d1, g1, r1,seed];
          disp(tmprun)
          res_cv = simures(n, d1, g1, r1,seed);
          
          inf(b,:) = [n, d1, g1, r1,seed];
          solve_time(b,:) = res_cv(1,:);
          run_time(b,:) = res_cv(2,:);
          RSS(b,:) = res_cv(3,:);
      
  end
  resn200 = [inf, solve_time, run_time, RSS];
  
  csvwrite("newsimudata/resn400_1.csv", [inf, solve_time, run_time, RSS])
  fprintf("n two is done")
  
  
  inf =  zeros(B,5);
  solve_time = zeros(B,3);
  run_time = zeros(B,3);
  RSS = zeros(B,3);

  parfor b = 1:B
              seed = 175342;
          seed = round(1.17*seed - b);
          n = 800;
          tmprun = [b,n, d1, g1, r1,seed];
          disp(tmprun)
          res_cv = simures(n, d1, g1, r1,seed);
          
          inf(b,:) = [n, d1, g1, r1,seed];
          solve_time(b,:) = res_cv(1,:);
          run_time(b,:) = res_cv(2,:);
          RSS(b,:) = res_cv(3,:);
      
  end
  resn200 = [inf, solve_time, run_time, RSS];
  
  csvwrite("newsimudata/resn800_1.csv", [inf, solve_time, run_time, RSS])
  fprintf("n three is done")
  
  inf =  zeros(B,5);
  solve_time = zeros(B,3);
  run_time = zeros(B,3);
  RSS = zeros(B,3);
  
  parfor b = 1:B 
              seed = 6632721;
          seed = round(0.243*seed - b);
          n = 1600;
          tmprun = [b,n, d1, g1, r1,seed];
          disp(tmprun)
          res_cv = simures(n, d1, g1, r1,seed);
          
          inf(b,:) = [n, d1, g1, r1,seed];
          solve_time(b,:) = res_cv(1,:);
          run_time(b,:) = res_cv(2,:);
          RSS(b,:) = res_cv(3,:);
      
  end
  resn2400 = [inf, solve_time, run_time, RSS];
  
  csvwrite("newsimudata/resn1600_1.csv", [inf, solve_time, run_time, RSS])
  
  fprintf("n is done")
  
  