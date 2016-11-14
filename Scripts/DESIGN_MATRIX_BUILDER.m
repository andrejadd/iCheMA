%
%
% Build the Design matrix given parent set, inhibition information and
% parameter vector 'K'.
%
%

function [D] = DESIGN_MATRIX_BUILDER(Data, parents, inhibition_vec, prod_ind, K)

    % The predictor data, i.e. the full design matrix containing all features
    data = Data.X;
    
    % number of data points.
    n_obs = size(data, 2);

    % The first covariate takes the degradation of the target node into account:
    D = (-1) * Data.X_degrad;

    
    % Add the data for each predictor to the design matrix:
    if (~isempty(parents))

        % Build additive covariates.
        if ( (prod_ind == 0) || (length(parents) == 1)) 

            % Generate a covariate vector for each parent node
            for i=1:length(parents)
                parent = parents(i);
                
                % Check if inhibitor or activitor
                if (inhibition_vec(i)==0) % if activator
                    new = data(parent,:)./(data(parent,:)+K(i)*ones(1,n_obs));   
                else % if inhibitor
                    new = (K(i)*ones(1,n_obs))./(data(parent,:)+K(i)*ones(1,n_obs));
                end

                D = [D;new];    
            end

        else % i.e. if prod_ind=1 AND length(parents)>1, then build a product as covariate

            % Generate ONE SINGLE PRODUCT covariate vector for all parent nodes
            new = ones(1,n_obs);

            for i=1:length(parents)
                parent = parents(i);
                if (inhibition_vec(i)==0)
                    new = new .* data(parent,:)./(data(parent,:)+K(i)*ones(1,n_obs));   
                else % activation
                    new = new .* (K(i)*ones(1,n_obs))./(data(parent,:)+K(i)*ones(1,n_obs));
                end
            end

            D = [D;new];
        end

    end

    D = D';

return
