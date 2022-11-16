function [lambda_seq] = lambda_gen_og_real(X,y,nlambda,opts)

lambda_max_init = 1e10;
beta1 = 0.9;
beta2 = 1.2;
%----------------------- Set optional items -----------------------
opts.tFlag=1;
 



flag  = 1;
while(flag)
    
   z = [0, lambda_max_init];
   [beta_est,~, ~, ~]= overlapping_LogisticR(X, y, z, opts);
   if(sum(abs(beta_est) > 1e-50) < 1e-20)
       lambda_max_init = lambda_max_init*beta1;
   else
       flag = 0;
   end
    
end

lambda_min_init =  1e-13;
flag  = 1;
while(flag)
    z = [0, lambda_min_init];
   [beta_est, c,~, ~]= overlapping_LogisticR(X, y, z, opts);

   
   if(sum(abs(beta_est) < 1e-50) < 1e-20)
       lambda_min_init = lambda_min_init*beta2;
   else
       flag = 0;
   end
end

C = (log(lambda_max_init) - log(lambda_min_init))/(nlambda - 1);
lambda_seq = lambda_min_init*exp(((1:nlambda) - 1)*C);
end
