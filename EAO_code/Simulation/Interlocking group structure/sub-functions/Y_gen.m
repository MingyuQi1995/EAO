%function to generate response variable
%Input: X (design matrix), beta (coefficient vector), sgm (variance of the error term), seed
%Output: a 1*n response vector.

function [Y] = Y_gen(X, beta, sgm,seed)
 
 rng(seed);
 n = size(X,1);
  
 Y  = X*beta +  normrnd(0,sgm,[1,n])';
      
end