%
% Copyright (c) 2016, Andrej Aderhold
%
% This script processes the results of the iCheMA main script and extract
% edge posterior probabilities. Given a known gold-standard (provided with
% the Biopepa example data) the AUROC and AUPREC scores are calculated
% with the edge probabilities.
%

function evaluate_Results()
    

    % This is for AUROC_AUPREC.m
    addpath('Scripts');

    % The results are located here 
    results_dir = sprintf('Results/');
    
    % The gradient type, can be 'RBFGradient' (analytic), or 'coarseGradient'.
    gradient_type = 'RBFGradient';  
    
    % A vector of data ids.
    data_instances = [1];

    % Different network names.
    biopepa_networks = {'wildtype'}; %, 'mutant_PRR7_PRR9', 'mutant_PRR5_PRR7_PRR9', 'mutant_PRR5_PRR7_PRR9_TOC1', 'mutant_TOC1', 'mutant_PRR7_PRR9_TOC1'};

    
    % Make AUROC and AUPREC calculation for each network and data instance
    % independently.
    for network = biopepa_networks

        network = biopepa_networks{1};
        
        % Loop over data instances
        for dataid = data_instances
            
            fprintf('handling network %s, data instance %i\n', network, dataid);

            edge_vec = process_Results(network, dataid, gradient_type, results_dir);
    
          
            % 
            % Get the gold-standard networks
            %
            if strcmp(network, 'wildtype')
                realnet = importdata('Data/Goldstandards/true_network_wildtype_mrna.txt');
            elseif strcmp(network, 'mutant_PRR7_PRR9')
                realnet = importdata('Data/Goldstandards/true_network_prr7-9_mrna.txt');
            elseif strcmp(network, 'mutant_PRR7_PRR9_TOC1')
                realnet = importdata('Data/Goldstandards/true_network_prr7-9-toc1_mrna.txt');
            elseif strcmp(network, 'mutant_PRR5_PRR7_PRR9')
                realnet = importdata('Data/Goldstandards/true_network_prr5-7-9_mrna.txt');
            elseif strcmp(network, 'mutant_PRR5_PRR7_PRR9_TOC1')
                realnet = importdata('Data/Goldstandards/true_network_prr5-7-9-toc1_mrna.txt');
            elseif strcmp(network, 'mutant_TOC1')
                realnet = importdata('Data/Goldstandards/true_network_toc1_mrna.txt');
            end

            % transpose realnet matrix and transform into vector, excluding -1 edges
            realnet = realnet';  % transpose, in order to make as.vector work and be consistent with weights.vec 
            realnet_vec = reshape(realnet, 1, size(realnet,1)*size(realnet,2));
            realnet_vec = realnet_vec(realnet_vec ~= -1);

            % Calculate the scores     
            [auroc, auprec] = AUROC_AUPREC(realnet_vec, edge_vec);
            
            % Save the edge posterior probabilities.
            fileout = sprintf(['%s/EVAL/EDGES_%s_%s_id%i.txt'], results_dir, network, gradient_type, dataid);
            edge_vec = edge_vec';
            dlmwrite(fileout, edge_vec, 'Delimiter', ',');

            % Save the AUROC and AUPREC scores
            fileout = sprintf(['%s/EVAL/AUROC_AUPREC_%s_%s_id%i.txt'], results_dir, network, gradient_type, dataid);
            auroc_auprec = [auroc, auprec];
            dlmwrite(fileout, auroc_auprec, 'Delimiter', ',');
        end
    end
end




%
% Extract the marginal log likelihoods for all the response and parent set
% combinations and calculate a posterior probability for each response
% parent pair. 
%
% Return the response - parent matrix as a vector. 
%
% 
function[edge_vec] = process_Results(network, dataid, gradient_type, results_dir)

    % Some basic settings to read in and handle the results data
    predictor_type = 'mRNA';
    RELOCATE_LIGHT_VAR = 1;  % Set this to 0 if not Biopepa data
    responses = 7;           % Biopepa specific setting   
    total_nodes = 8;         % Biopepa specific setting 
    fanin = 3;
    

    RESULTS = {};
    
    %
    %
    % Stage 1: Read in the marginal log likelihoods for each response and
    % each parent set.
    %
    %
    for response_node = 1:responses

        fprintf('.');

        % take out response_node as parent
        parent_nodes = setdiff(1:total_nodes, response_node);

        % this counter is used to identify the parent set generated with
        % recombination of different parent ids for different parent sizes
        parent_sets= {};
        parent_sets{end+1} = [];  % the first is the empty set

        max_parents = fanin;

        if fanin == -1
            max_parents = length(parent_nodes);
        end


        % create all possible parent set combinations from 1 parent node to
        % max_parents node.
        for ii = 1:max_parents

            comb_tmp = nchoosek(parent_nodes, ii);                                                                                                                     

            for j = 1:size(comb_tmp,1) 

                % append set
                parent_sets{end+1} = comb_tmp(j,:);

            end

        end         

        % Save it for later
        RESULTS{response_node}.parent_sets = parent_sets;
        RESULTS{response_node}.scores = {};


        response_max_val_logLL = -Inf;
        response_max_val_logLL_prior = -Inf;
       
        %
        % Loop over the different parent set configurations and read out
        % the marginal log likelihood (LL) for each one. Determine the highest
        % LL of each parent set given its different activator/inhibitor
        % setups. Ignore all other setups and use the highest one for
        % further probability calculations involving the parents in the
        % set.
        %
        for parent_set_id = 1:length(parent_sets)

            % Identify results file and load it
            filein = sprintf('%s/%s/OUT_r%i_psi%i_BIOPEPA-%s_run%i.mat', results_dir, gradient_type, response_node, parent_set_id, network, dataid);
            load(filein);

            % Save the marginal log likelihood
            RESULTS{response_node}.scores{parent_set_id}.scores = results;

            % Use this to identify the parent set with the highest log
            % likelihood
            score_tmp = results;

            % Just to be save, delete results.
            clear results;

            %
            % get max values for each parent_set
            %
            max_val_logLL = -Inf;

            % loop through cell array and extract largest likelihood
            for jj = 1:length(score_tmp)

                if score_tmp{jj} > max_val_logLL
                    max_val_logLL = score_tmp{jj};
                end
            end

            RESULTS{response_node}.scores{parent_set_id}.max_score = max_val_logLL;


            % Check if this is the overall maximum score for this response
            % (used to move likelihood to maximum of 0 for better exp()
            % translation
            if max_val_logLL > response_max_val_logLL
                response_max_val_logLL = max_val_logLL;
            end


        end

        RESULTS{response_node}.response_max_score =  response_max_val_logLL; 

    end  % End of response

    fprintf('\n');

    
    
    %
    %
    % Stage 2: For each response, calculate the probability of the parents.
    %
    %

    SCORES = NaN * zeros(responses, total_nodes);
    SCORES_ALL = NaN * zeros(responses, total_nodes);

    for response_node = 1:responses

        result = RESULTS{response_node};

        % for each possible parent
        for parent = 1:total_nodes

            % ignore self-loop
            if parent == response_node
                continue;
            end

            % these are used to calculate probability, marginalized over all
            % scores (score_edge_all)
            score_edge_logLL = 0;
            score_edge_all_logLL = 0;

            % look into all parent_sets
            for parent_set_id = 1:length(result.parent_sets)

                % --> max value for parent set (can have different
                % values dependent on activator/inhibitor setup
                % result.scores{parent_set_id}.max_Thermo_logLL_prior
                % 
                % --> overall max. value for all parent set, used to
                % shift all log score to max value of 0 for better
                % exp() usage
                % result.response_max_Thermo_logLL_prior

                % only likelihood
                tmp_ll_logLL = result.scores{parent_set_id}.max_score - result.response_max_score;
                tmp_prob_logLL = exp(tmp_ll_logLL);


                % if the parent is in the parent_set .. take the intercept of
                % parent set and the parent, if it is 0 then its in the set,
                % negate it and we are in side the if condition
                if ~isempty(intersect(result.parent_sets{parent_set_id}, parent))

                    score_edge_logLL = score_edge_logLL + tmp_prob_logLL;

                end

                score_edge_all_logLL = score_edge_all_logLL + tmp_prob_logLL;

            end

            SCORES(response_node, parent) = score_edge_logLL;
            SCORES_ALL(response_node, parent) = score_edge_all_logLL;

        end
    end

    
    DAG_PROB = SCORES ./ SCORES_ALL;
    DAG_PROB(isnan(DAG_PROB)) = -1;

    %
    % Relocate light*P to 3rd position.
    % 
    if(RELOCATE_LIGHT_VAR)
        DAG_PROB = [DAG_PROB(:,1:2),DAG_PROB(:,8), DAG_PROB(:,3:7)];
    end
    
    
    % Transpose and convert to vector.
    DAG_PROB = DAG_PROB';
    edge_vec = DAG_PROB(:);

    % Get rid of -1 (self-loops). This variable is returned.
    edge_vec( find(edge_vec == -1) ) = [];

end

    