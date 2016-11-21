
%
% Given n available numbers and r=(0,..,max_r) possible combinations of these numbers at a time (order does not matter), 
% calculate the total amount of possible combinations including the empty set. 
 %      
function [n_comb] = n_combinations(n, max_r)

    n_comb = 0;

    % precompute
    n_fact = factorial(n);
    
    %
    % sum(C(n,r)) , i.e. add up all the combinations for n available
    % numbers and r combination set sizes. 
    %
    for r = 0:max_r
        tmp = n_fact / (factorial(n-r)*factorial(r));
        n_comb = n_comb + tmp;
    end
    
end