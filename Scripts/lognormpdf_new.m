
function [log_densities] = lognormpdf_new(vector,mu,std_dev)

n = length(vector);

log_densities = zeros(1,n);

for i=1:n
    
    value = vector(i);
    
    log_densities(1,i) = -0.5 * log(2*pi) - log(std_dev) - 0.5/(std_dev^2) * (value-mu)^2; 
        
end

return
