%Function to generate groups 
%Input: d (group size), g (number of groups) and r (overlapping ratio between two groups)
%Output: A g*p matrix where each row represents a specific group.


function [G] = Gint_gen(d, g, r)

overlap_num = ceil(d*r);

p0 = d*g;

G = zeros(g, p0);

G(1,1:d) =  1;

k = d - overlap_num + 1;

 for i = 2:g
       k1 = k + d - 1;
    G(i, k:k1) = 1;
    k = k1 - overlap_num + 1 ;
    end

 G = G(:,1:k1);
  
end