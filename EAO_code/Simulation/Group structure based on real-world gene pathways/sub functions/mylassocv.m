function [Delta, lambda_min] = mylassocv(X, y, lambda_seq, opts, beta_real)



nlambda = length(lambda_seq);

Delta = [];
%st = [];

opts.init=2;   

rho = lambda_seq(nlambda);

%tic;
[beta_est, ~, ~]= LeastR(X, y, rho, opts);
%stime = toc;
%st(nlambda) =stime;

Delta(nlambda) =  sum((beta_est - beta_real).^2)./(sum(beta_real.^2));

 
    for j = 1:(nlambda-1)
        opts.init = 1;
        opts.x0 = beta_est;
        rho = lambda_seq(nlambda - j);
       
       % tic;
        [beta_est, ~, ~]= LeastR(X, y, rho, opts);
       % stime = toc;
       % st(nlambda - j) = stime;
     %   y_est = X_test*beta_est;
        
%        cv_res(i,j) = sum((y_est - y_test).^2);
       Delta(nlambda - j) =  sum((beta_est - beta_real).^2)./(sum(beta_real.^2));
 
    end

Delta_min = min(Delta);
idx_min = find(Delta == Delta_min,1,"last");
lambda_min = lambda_seq(idx_min);
%end
%cv_m = mean(cv_res,1);
%cv_min = min(cv_m);
%idx_min = find(cv_m == cv_min);
%lambda_min = lambda_seq(idx_min);
end