%Function of getting the index to apply SLEP
%Input: non-overlapping groups
%Output: Index for each group


function [ind] = indgroup(group)
  ind = [];
  k = length(group) - 1;
  
  
 for i = 1:k
     if(group(i) ~= group(i+1))
         ind = [ind,i] ;
     end
 end
  ind = [0,ind,length(group)];
  
end