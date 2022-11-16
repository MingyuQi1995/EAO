function [lambda_seq] = reallassolambda_gen(X,y,nlambda,opts)

lambda_max_init = 1e8;
beta1 = 0.9;
beta2 = 1.2;
%----------------------- Set optional items -----------------------
opts.tFlag=1;


flag  = 1;
opts.tFlag=0;
while(flag)
    
   rho = lambda_max_init;
   [beta_est,~, ~, ~]=  LogisticR(X, y, rho, opts);
   if(sum(abs(beta_est) > 1e-50) < 1e-20)
       lambda_max_init = lambda_max_init*beta1;
   else
       flag = 0;
   end
    
end

lambda_min_init = 1e-12;
flag  = 1;
while(flag)
    rho =  lambda_min_init;
   [beta_est,~, ~, ~]= LogisticR(X, y, rho, opts);
   if(sum(abs(beta_est) < 1e-50) < 1e-20)
       lambda_min_init = lambda_min_init*beta2;
   else
       flag = 0;
   end
end

C = (log(lambda_max_init) - log(lambda_min_init))/(nlambda - 1);
lambda_seq = lambda_min_init*exp(((1:nlambda) - 1)*C);
end