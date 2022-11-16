%Function to use a line search to find the regularization path for non-overlapping group lasso
%Additionally, this function generates on a log scale similar to glm package in R.
%Input: X,y, nlambda (the number of elements in the path), opts (see SLEP: Sparse Learning with Efficient Projections)
%Output: A sequence of lambda 

function [lambda_seq] = intlambda_gen(X,y,nlambda,opts)

lambda_max_init = 1e8;
beta1 = 0.8;
beta2 = 1.2;
%----------------------- Set optional items -----------------------



flag  = 1;
while(flag)
    
   rho = lambda_max_init;
   [beta_est, ~, ~]= glLeastR(X, y, rho, opts);
   if(sum(abs(beta_est) > 1e-5) < 1e-9)
       lambda_max_init = lambda_max_init*beta1;
   else
       flag = 0;
   end
    
end

lambda_min_init = 1e-8;
flag  = 1;
while(flag)
    rho =  lambda_min_init;
   [beta_est, ~, ~]= glLeastR(X, y, rho, opts);
   if(sum(abs(beta_est) < 1e-5) < 1e-9)
       lambda_min_init = lambda_min_init*beta2;
   else
       flag = 0;
   end
end

C = (log(lambda_max_init) - log(lambda_min_init))/(nlambda - 1);
lambda_seq = lambda_min_init*exp(((1:nlambda) - 1)*C);
end