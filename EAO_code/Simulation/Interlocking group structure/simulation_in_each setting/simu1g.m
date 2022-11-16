   B = 100;
  
  n1 = 800;
  r1=0.2;
  d1=20;
  inf =  zeros(B,5);
  solve_time = zeros(B,3);
  run_time = zeros(B,3);
  RSS = zeros(B,3);
    parpool('local',100)
  parfor b = 1:B
            seed = 57333;
          seed = round(9.303*seed - b);
         g = 200;
           tmprun = [b,n1, d1, g, r1,seed];
           disp(tmprun)
         res_cv = simures(n1, d1, g, r1,seed);
          inf(b,:) = [n1, d1, g, r1,seed];
          solve_time(b,:) = res_cv(1,:);
          run_time(b,:) = res_cv(2,:);
          RSS(b,:) = res_cv(3,:);
     
  end
  resg100 = [inf,solve_time, run_time, RSS];
   csvwrite("newsimudata/resg200_1.csv", [inf, solve_time, run_time, RSS])
   fprintf("g one is done")
   
    
  

  inf =  zeros(B,5);
  solve_time = zeros(B,3);
  run_time = zeros(B,3);
  RSS = zeros(B,3);

  parfor b = 1:B
            seed =  17918;
          seed = round(20.22*seed - b);
         g = 400;
           tmprun = [b,n1, d1, g, r1,seed];
           disp(tmprun)
         res_cv = simures(n1, d1, g, r1,seed);
          inf(b,:) = [n1, d1, g, r1,seed];
          solve_time(b,:) = res_cv(1,:);
          run_time(b,:) = res_cv(2,:);
          RSS(b,:) = res_cv(3,:);
      
  end
  resg200 = [inf,solve_time, run_time, RSS];
   csvwrite("newsimudata/resg400_1.csv", [inf, solve_time, run_time, RSS])
   fprintf("g two is done")
   
   
   
  inf =  zeros(B,5);
  solve_time = zeros(B,3);
  run_time = zeros(B,3);
  RSS = zeros(B,3);
 
  parfor b = 1:B
            seed = 228188;
          seed = round(0.29*seed - b);
         g = 800;
           tmprun = [b,n1, d1, g, r1,seed];
           disp(tmprun)
         res_cv = simures(n1, d1, g, r1,seed);
          inf(b,:) = [n1, d1, g, r1,seed];
          solve_time(b,:) = res_cv(1,:);
          run_time(b,:) = res_cv(2,:);
          RSS(b,:) = res_cv(3,:);
      
  end
  resg200 = [inf,solve_time, run_time, RSS];
   csvwrite("newsimudata/resg800_1.csv", [inf, solve_time, run_time, RSS])
   fprintf("g two is done")

   

  inf =  zeros(B,5);
  solve_time = zeros(B,3);
  run_time = zeros(B,3);
  RSS = zeros(B,3);

  parfor b = 1:B 
            seed = 9971482;
          seed = round(5.142*seed - b);
         g = 100;
           tmprun = [b,n1, d1, g, r1,seed];
           disp(tmprun)
         res_cv = simures(n1, d1, g, r1,seed);
          inf(b,:) = [n1, d1, g, r1,seed];
          solve_time(b,:) = res_cv(1,:);
          run_time(b,:) = res_cv(2,:);
          RSS(b,:) = res_cv(3,:);
     
  end
  resg800 = [inf,solve_time, run_time, RSS];
   csvwrite("newsimudata/resg100_1.csv", [inf, solve_time, run_time, RSS])
   fprintf("g  is done")
   