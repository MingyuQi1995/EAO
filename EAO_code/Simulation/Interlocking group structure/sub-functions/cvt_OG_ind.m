% function to get the group information to apply SLEP package
% Input: G (group structure)
% Output: opt_G and opt_ind (see SLEP: Sparse Learning with Efficient Projections for details)

function [opt_G,opt_ind] = cvt_OG_ind(G)
  p =  size(G,2);
  g =  size(G,1);
  opt_G = [];
  
 opt_ind = zeros(3, g);

  l0 = 1;
  
 for i = 1:g
     tmp_idx = find(G(i,:) == 1);
  
    l1 = l0 + length(tmp_idx) ;
    opt_G = [opt_G, tmp_idx];
    opt_ind(1,i) = l0;
    opt_ind(2,i) = l1-1;
    opt_ind(3,i) = 1;
    l0 = l1;
 end
  
end