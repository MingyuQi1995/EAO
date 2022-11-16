%Function to partition overlapping groups into non-overlapping groups
%Input: G (group structure)
%Output: Group (non-overlapping groups), w (weight for each non-overlapping group)


function [Group,w] = gen_par(G)
  G1 = mypar(G);
  H = sum(G);
  k = size(G1,1);
  col =  size(G1,2);
  H1 = [];
  
 for i = 1:k
     tmp_idx = find(G1(i,:) == 1);
    H1(i)  = mean(H(tmp_idx));

 end
   w = H1;

  Group = [];
  
  for i = 1:col
     idx = find(G1(:,i) == 1);
    Group(i)   = idx;

 end
  
  
end