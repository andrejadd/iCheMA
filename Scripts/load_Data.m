
%
% load data
%
    
function[X, y, X_degrad, light_var] = load_Data(runid, response_node, ...
                                         biopepa_network, gradient, predictor_type)

        

    %
    % load mRNA data  -- note, we need to load both, i.e. the mRNA
    % for the self-loop is also required for the protein predictors 
    %
    loadfile = sprintf(['Data/biopepa_%s_%s_mRNA-predictors_id%i.csv'], ...
                       biopepa_network, gradient, runid);
    X_mRNA = importdata(loadfile);
    X_mRNA = X_mRNA.data;
    
    % put light*P (at column 3) as last entry, to make it easier for the parent-set
    % generation process

    X_mRNA = [X_mRNA, X_mRNA(:,3)]; % for mRNAs
    light_var = X_mRNA(:,3);  % save for seperate use
    X_mRNA(:,3) = [];

    % need this always, especially when we have proteins as predictors
    X_degrad = X_mRNA(:,response_node)';
    
    
    % 
    % load protein data
    %
    loadfile = sprintf(['Data/biopepa_%s_%s_protein-predictors_id%i.csv'], ...
                       biopepa_network, gradient, runid);
        
    X_protein = importdata(loadfile);
    X_protein = X_protein.data;
    
    % put light*P (at column 4) as last entry, to make it easier for the parent-set
    % generation process
    X_protein = [X_protein, X_protein(:,4)]; % for proteins
    X_protein(:,4) = [];

    
    %
    % pick the one that is required
    %
    if strcmp(predictor_type, 'protein')
        X = X_protein;
    else
        X = X_mRNA;
    end
    
    
    
    % load the response (gradient) data
    loadfile = sprintf(['Data/biopepa_%s_%s_gradient-' ...
                        'response_id%i.csv'], biopepa_network, gradient, runid);
    Y = importdata(loadfile);
    Y = Y.data;
    y = Y(:,response_node);
 
    

    %
    % Filter out data where target is knocked-out  
    %
    %  everytime GI, LHY, TOC1, PRR7 or PRR9 are a response (target) it is necessary to exlude parts of the observations that
    %            have been knock-out for the corresponding gene, e.g.
    %
    %   Target  ID  Positions
    %   GI       2    1:13
    %   LHY      4   14:26
    %   TOC1    14   27:39
    %   PRR7    10   40:52
    %   PRR9    12   40:52
    %
    %   The PRR's are in a double mutant, so they share the same Position in the design matrix
    %
    %   Do this: Tf a target/response variable is encountered that matches 'ID', exclude corresponding 'Positions' from all calculations !
    %            Only do this for the constant knock-outs but not the pruned networks, which we assume to not know
    %

    exclude_obs = [];
    
    % take out observations if the current response is one of the below and if this feature is enabled
    if(response_node == 1) exclude_obs = 1:13;  end
    if(response_node == 2) exclude_obs = 14:26; end
    if(response_node == 6) exclude_obs = 27:39; end
    if(response_node == 4) exclude_obs = 40:52; end
    if(response_node == 5) exclude_obs = 40:52; end
        
    nr_obs = size(X,1);
    
    take_these_obs = setdiff(1:nr_obs, exclude_obs);

    % take only the ones where the target is not ko'd
    X = X(take_these_obs,:);
    X = X';

    X_degrad = X_degrad(take_these_obs);
    
    y = y(take_these_obs);

        

end

