%Function to generate the design matrix
%Input: n (sample size), d (group size), G (a g*p group matrix), toep_cor (the correlation within the same group), seed (arbitrary value to reproduce the result)
%Output: a n*p design matrix.

function [X] = data_gen(n, d, G,toep_cor, seed)
    
    rng(seed);
    
    p = size(G,2);
    g = size(G,1);
    
    mu =   zeros(1,p);
    toep = zeros(p,p);
    G1 = mypar(G);
    k = size(G1,1);
     for i = 1:g
     idx_tmp = find(G(i,:) == 1);
    toep(idx_tmp, idx_tmp) = 0.36;
     end
 
      for i = 1:k
     idx_tmp = find(G1(i,:) == 1);
    toep(idx_tmp, idx_tmp) = toep_cor;
      end
     
   Sigma = toep;  
   Sigma = Sigma + (1 - toep_cor ) *diag(repelem(1 , p));
   [V, D] = eig(Sigma);
   D1 = diag(D);
   D1(D1<0.1) = 0.1;
   D = diag(D1);
   Sigma = V*D*V';
  

   X = mvnrnd(mu,Sigma,n);
      
end