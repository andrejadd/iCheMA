% Copyright 2016, Marco Grzegorczyk and Andrej Aderhold
%
% Function run_ICHEMA()
%
% This is the Matlab implementation of iCheMA with adaptations to the Biopepa data
% set. See the README.txt for more information.
%
% Please cite the following paper when using this software:
%
%  Aderhold, A., Grzegorczyk, M., and Husmeier, D. (2016). Approximate
%  Bayesian inference in semi-mechanistic models. Statistics and Computing, 1-38.
%
%

function run_ICHEMA()


    %
    % Example of using this iCheMA implementation with the Biopepa data.
    %   This data has 7 response nodes. Assuming a max_fanin of 3 this
    %   results in 65 possible parent set configuration for each response.
    %
    %  Here we loop over each response and each of the possible parent set
    %  configuration. See the function getParentSets() to understand how
    %  the parent sets are generated. This should be adapted your own
    %  needs.
    %
    % Note: Running this with a sufficiently high number of iterations
    % can take a long time because of all the possible parent set configurations 
    % and activation and inhibition term types. Consider Running this on a cluster
    % is ...
     
    MCMC_iterations = 300;
    
    max_response_nodes = 7;  % Biopepa specific. 
    max_parent_set_ids = 64; % Given 7 response nodes, possible number of parent configurations.
    
    for response_node = 1:max_response_nodes
        
        for parent_set_id = 1:max_parent_set_ids
    
            % Run iCheMA, at the end of this function the results are saved
            % into the 'Results/' directory.
            %
            % See the following parameter settings inside run_ICHEMA_core() that define the data to be read:
            %    gradient 
            %    predictor_type
            %    network
            %    data_instance  
            %
            run_ICHEMA_core(response_node, parent_set_id, MCMC_iterations)
        
        end
        
    end
end

%
% Main entry method: This function reads the data, generates the parent
% set, the activation/inhibition flag vector, and finally calls  
%
%
%
function[] = run_ICHEMA_core(response_node, parent_set_id, iterations)

    % All the script are in here:
    addpath('Scripts');

    rng('shuffle');

    % Define the input data: This information is used in load_Data() to
    % identify the predictor data (network, predictor type, data instance)
    % and the response data (gradient, network, data instance)
    %
    % NOTE: Adapt these and the content of load_Data() for you own needs.
    gradient = 'RBFGradient';  % the response type
                               % 'RBFGRadient' : Analytic gradient
                               % 'coarseGradient' : Numerical gradient
    predictor_type = 'mRNA';   % identify the predictor type
    network = 'wildtype';      % identify the network
    data_instance = 1;         % one of the five data instances
    
    % Maximum number of child nodes, default setting
    % When changing this the parent set function also becomes affected, so
    % threat with care and adapt getParentSets() accordingly.
    max_fanin = 3;
   
    % The results are written into here
    fileout = sprintf('Results/OUT_r%i_psi%i_BIOPEPA-%s_run%i.mat', response_node, parent_set_id, network, data_instance);

    % This function loads the design matrix 'X', the response vector 'y'
    % (depends on variable 'response_node'), the degradation term
    % 'X_degrad', and the light indicator vector 'light_var'. The former vector
    % is specific to the data used here.
    %
    % IMPORTANT: Adapt this function to your own requirements. This is an
    % example for the Biopepa data as it was used in the paper.
    [X, y, X_degrad, light_var] = load_Data(data_instance, response_node, network, gradient, predictor_type);

    % Pack it together into one data structure.
    Data.X = X;
    Data.X_degrad = X_degrad;
    Data.y = y;
    Data.light_var = light_var;
    
    % Number of nodes
    total_nodes = size(Data.X, 1);

    % Take out response_node as parent, so this is a vector with putative parent nodes.
    parent_nodes = setdiff(1:total_nodes, response_node);
           
    % Map the parent_set_id to a vector of parent nodes and a corresponding
    % matrix with activation/inhibition combinations. 
    %
    % E.g. If parent_set_id is 10, the number of total_nodes is 8, the 
    %      response_node is 2, and max_fanin is 3, then parent_nodes = (1,3,4,5,6,7,8) 
    %      and parent_set is [1,4] with inhib_mat = [[0,0];[0,1];[1,0];[1,1]]
    %
    %   This means that calculations for response 2 will be done, assuming
    %   the parent nodes are 1 and 4, and that there are four different
    %   combination of activation and inhibition.
    %   
    %  The first row in inhib_mat is [0,0], meaning that 1 and 4 are
    %  both activator. The second row is [0,1], meaning 1 is an activator
    %  and 4 is an inhibitor, etc..
    %
    % NOTE: Modify this function to you own needs:
    [parent_set, inhib_mat] = getParentSets(parent_nodes, parent_set_id, max_fanin);
       
    % Put the scores here - This is saved to the output file.
    results = {};
    
    % If the parent set is empty, call the main function differently that
    % if the parent set is not empty. 
    if isempty(parent_set)

        % Call functions with lack of any parent nodes, i.e. parent_set is
        % an empty vector, and the argument that follows is also empty ('[]'). 
        fprintf('Compute marginal log likelihood with iCheMa for empty parent set..\n');

        results{end+1} = ICHEMA_marginalLL(Data, response_node, parent_set, [], iterations);
       
    else  % .. parent set is not empty

        % Loop over the different activiation/inhibition combinations
        % stored in the matrix 'inhib_mat'
        for ii = 1:size(inhib_mat, 1)

            % get the activation/inhibition flag vector
            inhib_set = inhib_mat(ii,:);
            fprintf('Compute marginal log likelihood with iCheMa for activation/inhibition flag set %i..\n', ii);
            results{end+1}  = ICHEMA_marginalLL(Data, response_node, parent_set, inhib_set, iterations);
                
        end
    end

    fprintf('\nattempting to write to %s\n', fileout);
    save(fileout, 'results', 'parent_set', 'inhib_mat', 'parent_set_id', 'response_node', 'max_fanin');   
    

end

