%Function to generate coefficient vector
%Input: G (group structure), seed
%Output: a 1*p coefficient vector

function [beta] =beta_gen0(G, seed)
    rng(seed);
    p = size(G,2);
    g = size(G,1);
    beta  = normrnd(10,4,[1,p]);
    randomNumbers = randi([0, 1], [1, p]);
    ran = (randomNumbers - 0.5)*2;
    beta = beta .* ran;
    
    k = floor(0.91*g);
     
    zerogroup = randsample([1:g],k,false);
    
 
     for i = 1:length(zerogroup)
    idx_tmp = find(G(zerogroup(i),:) == 1);
    beta(idx_tmp) = 0;
     end
 
 
      
end