
B = 100;



 
  
  
  
   r1= 0.2;
   g1 = 200;
  n1=  800;
  inf =  zeros(B,5);
  solve_time = zeros(B,3);
  run_time = zeros(B,3);
  RSS = zeros(B,3);
    parpool('local',33)
  parfor b = 1:B
          seed = 24952992;
          seed = round(1.757*seed - b);
         d = 10;
           tmprun = [b,n1, d, g1, r1,seed];
            disp(tmprun)
          res_cv = simures(n1, d, g1, r1,seed);
          
          
          
          inf(b,:) = [n1, d, g1, r1,seed];
          solve_time(b,:) = res_cv(1,:);
          run_time(b,:) = res_cv(2,:);
          RSS(b,:) = res_cv(3,:);
      
  end
  
    csvwrite("newsimudata/resd10_1.csv", [inf, solve_time, run_time, RSS])
   fprintf("d one is done")
   
  inf =  zeros(B,5);
  solve_time = zeros(B,3);
  run_time = zeros(B,3);
  RSS = zeros(B,3);
 
  parfor b = 1:B
             seed = 3909599;
          seed = round(15.625*seed - b);
          d = 20;
           tmprun = [b,n1, d, g1, r1,seed];
          disp(tmprun)
          res_cv = simures(n1, d, g1, r1,seed);
          
          
          
          inf(b,:) = [n1, d, g1, r1,seed];
          solve_time(b,:) = res_cv(1,:);
          run_time(b,:) = res_cv(2,:);
          RSS(b,:) = res_cv(3,:);
      
  end
    csvwrite("newsimudata/resd20_1.csv", [inf, solve_time, run_time, RSS])
   fprintf("d two is done") 



  inf =  zeros(B,5);
  solve_time = zeros(B,3);
  run_time = zeros(B,3);
  RSS = zeros(B,3);
  parfor b = 1:B
        seed =  205902;
          seed = round(6.171*seed - b);
          d = 40
           tmprun = [b,n1, d, g1, r1,seed];
        disp(tmprun)
          res_cv = simures(n1, d, g1, r1,seed);
          
          inf(b,:) = [n1, d, g1, r1,seed];
          solve_time(b,:) = res_cv(1,:);
          run_time(b,:) = res_cv(2,:);
          RSS(b,:) = res_cv(3,:);
      
  end
  
   csvwrite("newsimudata/resd40_1.csv", [inf, solve_time, run_time, RSS])
   fprintf("d three is done")

  inf =  zeros(B,5);
  solve_time = zeros(B,3);
  run_time = zeros(B,3);
  RSS = zeros(B,3);

  parfor b = 1:B
             seed = 112342;
          seed = round(12.134*seed - b);
          d = 80;
           tmprun = [b,n1, d, g1, r1,seed];
          disp(tmprun)
          res_cv = simures(n1, d, g1, r1,seed);
          
          
          
          inf(b,:) = [n1, d, g1, r1,seed];
          solve_time(b,:) = res_cv(1,:);
          run_time(b,:) = res_cv(2,:);
          RSS(b,:) = res_cv(3,:);
      
  end
  

  
   csvwrite("newsimudata/resd80_1.csv", [inf, solve_time, run_time, RSS])
   fprintf("d  is done")
   
   