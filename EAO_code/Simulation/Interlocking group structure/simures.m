%Main function for simulation
%Input: n (sample size), d (group size), g (number of groups), r (overlapping ratio), seed (an arbitrary number for reproducing same results)
%Output: Res (Regularization path computing time and estimation error of each method under the setting).


function[Res] = simures(n, d, g, r, seed)
 
  G = Gint_gen(d, g, r); 
  
  X =  data_gen(n, d, G, 0.6, seed);
  
  beta_real = beta_gen0(G,seed)';
  
  Expect = X*beta_real;
  
  sig = mean(abs(Expect));
  
  noi =   sig/3;

  y = Y_gen(X, beta_real,  noi, seed);  


  root = pwd;
  root1 = strcat(root, "/SLEP-master");
  addpath(genpath("./SLEP-master"));


%----------------------- Set optional items -----------------------
opts=[];

% Starting point
opts.init=2;        % starting from a zero point

% Termination 

opts.tFlag=0;
opts.maxIter=3e+5;  % maximum number of iterations
opts.tol= 1e-3;      % the tolerance parameter

% regularization
opts.rFlag=0;       % use ratio

% Normalization
opts.nFlag=0;       % without normalization

[opt_G,opt_ind] = cvt_OG_ind(G);
opts.G = opt_G;
opts.ind = opt_ind;

opts.maxIter2 = 3e+5;
opts.tol2= 1e-4;
opts.flag2=2;

% regularization
opts.rFlag=0;       % the input parameter 'rho' is a ratio in (0, 1)

opts.rsL2=0;

opts.mFlag=1;       % treating it as compositive function 
opts.lFlag=1;       % Nemirovski's line search

nlambda = 100;

tic;
lambda_seq = soglambda_gen(X,y,nlambda,opts);
[cv_res,lambda_min] = myoglcv(X, y, lambda_seq, opts, beta_real);
solve_time_og = toc;
solve_time_og = solve_time_og;

opts.init=2;  

z = [0, lambda_min];

tic;
[beta_est_og, ~, ~]=  overlapping_LeastR(X, y, z, opts);
run_Time_og = toc;

RSSog = sqrt(sum((beta_est_og - beta_real ).^2))./sqrt((sum(beta_real.^2)));
OG_res = [solve_time_og;run_Time_og;RSSog;beta_est_og];

  
%----------------------- wlas -----------------------






lasweight = sum(G)';
beta_realwlas = beta_real.* lasweight;
Xwlas = X ./(lasweight');

nlambda = 100;


tic;

lambda_seq = slassolambda_gen(Xwlas,y,nlambda,opts);

[cv_res,lambda_min] = mylassocv(Xwlas, y, lambda_seq, opts, beta_realwlas);
solve_time_wlas = toc;

opts.init=2;  

rho = lambda_min;

tic;
[beta_est_wlas, ~, ~]=  LeastR(Xwlas, y, rho, opts);
run_Time_wlas = toc;

beta_est_wlas = beta_est_wlas./lasweight;

RSS_wlas = sqrt(sum((beta_est_wlas - beta_real).^2))./sqrt((sum(beta_real.^2)));

wlas_res = [solve_time_wlas;run_Time_wlas;RSS_wlas;beta_est_wlas];



%% wglasso

[Group,w] = gen_par(G);
ind =  indgroup(Group);    % number of groups
q=2;                 % the value of q in the L1/Lq regularization


% Group Property


opts.ind=ind.';       % set the group indices
opts.q=q;           % set the value for q % set the weight for positive and negative samples
opts.gWeight=w'; % set the weight for the group, a cloumn vector


opts.mFlag=0;       % treating it as compositive function 
opts.lFlag=0;       % Nemirovski's line search


nlambda = 100;

tic;

lambda_seq = intlambda_gen(X,y,nlambda,opts);
[cv_res,lambda_min] = intmywglcv(X, y, lambda_seq, opts, beta_real);
solve_time_wgl = toc;
%lambda_min = lambda * n;

opts.init=2;  
z = lambda_min;


tic;
[beta_est_wgl, val1, val2]= glLeastR(X, y, z, opts);
run_Time_wgl = toc;


RSS_wgl = sqrt(sum((beta_est_wgl -  beta_real).^2))./sqrt((sum( beta_real.^2)));
wgl_res = [solve_time_wgl;run_Time_wgl;RSS_wgl;beta_est_wgl];


%% Final
Res = [OG_res wlas_res wgl_res];
end
  