function[canres] = simucan(seed)
%%  Oglasso


 root = pwd;
  root1 = strcat(root, "/SLEP-master");
  addpath(genpath("./SLEP-master"));

load('realsimudata/simudatacan.mat');

X = data.X;
opt_G = data.opt_G;
opt_ind = data.opt_ind;

rng(seed);

    p = size(X,2);
    beta  = normrnd(10,4,[1,p]);
    randomNumbers = randi([0, 1], [1, p]);
    ran = (randomNumbers - 0.5)*2;
    beta = beta .* ran;
    
    ksize = length(opt_ind(3,:));
     kks = 30;
   indtemp = find(opt_ind(3,:) < 20);
   Smpl = randsample(indtemp, kks);
   optzero = opt_ind;
   optzero(:,Smpl) = [];
   for i = 1:length(optzero)
    idx_tmpa = optzero(1,i);
    idx_tmpb = optzero(2,i);
    idx_tmp = opt_G(idx_tmpa:idx_tmpb);
    beta(idx_tmp) = 0;
   end
   
   beta_real = beta';

     Expect = X*beta_real;
  
  sig = mean(abs(Expect));
  
  noi =   sig/3;
  
  noisd = sqrt(noi);

  y = Y_gen(X, beta_real,  noisd,seed);  
  
  
  opts=[];

% Starting point
opts.init=2;        % starting from a zero point

% Termination 

opts.tFlag=0;
opts.maxIter=3e+5;  % maximum number of iterations
opts.tol= 1e-4;      % the tolerance parameter

% regularization
opts.rFlag=0;       % use ratio

% Normalization
opts.nFlag=0;       % without normalization

opts.G = opt_G;
opts.ind = opt_ind;

opts.maxIter2 = 3e+5;
opts.tol2= 1e-4;
opts.flag2=2;

% regularization
opts.rFlag=0;       % the input parameter 'rho' is a ratio in (0, 1)

opts.rsL2=0;

opts.mFlag=0;       % treating it as compositive function 
opts.lFlag=0;       % Nemirovski's line search

nlambda = 100;

tic;
lambda_seq = oglambda_gen(X,y,nlambda,opts);
[cv_res,lambda_min] = myoglcv(X, y, lambda_seq, opts, beta_real);
solve_time_og = toc;


opts.init=2;  

z = [0, lambda_min];

tic;
[beta_est_og, ~, ~]=  overlapping_LeastR(X, y, z, opts);
run_Time_og = toc;

RSSog = sqrt(sum((beta_est_og - beta_real ).^2))./sqrt((sum(beta_real.^2)));
OG_res = [solve_time_og;run_Time_og;RSSog;beta_est_og];
  
%% wlasso
load('realsimudata/simulassodatacan.mat');

lasweight = data.W;
beta_realwlas = beta_real.* lasweight;
Xwlas = X ./(lasweight');

nlambda = 100;


lambda_seq = lassolambda_gen(Xwlas,y,nlambda,opts);

tic;
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


%%  wsglasso

load('realsimudata/simuwsgdatacan.mat');

X = data.X;

opts.ind = data.opt_ind;



tic;

lambda_seq = lambda_gen(X,y,nlambda,opts);
[cv_res,lambda_min] = mywglcv(X, y, lambda_seq, opts, beta_real);
solve_time_wsg = toc;

%lambda_min = lambda * n;
opts.init=2;  
z = [0, lambda_min];


tic;
[beta_est_wsg, ~, ~]= sgLeastR(X, y, z, opts);
run_Time_wsg = toc;

RSS_wsg = sqrt(sum((beta_est_wsg -  beta_real).^2))./sqrt((sum( beta_real.^2)));
wsg_res = [solve_time_wsg;run_Time_wsg;RSS_wsg;beta_est_wsg];

canres = [OG_res wlas_res wsg_res];
end