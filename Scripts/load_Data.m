
%
% Load the data. This includes loading the predictor data and the response
% data. This function also does some Biopepa specific data handling if the
% function argument BIOPEPA_DATA is set to 1 (true). 
%
% Returns: 
%
%   Data.X         - Design matrix with all predictor variables
%   Data.X_degrad  - The concentration data for the response node (self-loop)
%   Data.y         - response data as change of concentration (gradient)
%   Data.light_vec - a vector with light information, Biopepa specific
%
    
function[Data] = load_Data(runid, response_node, biopepa_network, gradient, BIOPEPA_DATA)

        

    %
    % load mRNA data  -- note, we need to load both, i.e. the mRNA
    % for the self-loop is also required for the protein predictors 
    %
    loadfile = sprintf(['Data/biopepa_%s_mRNA-predictors_id%i.csv'], ...
                       biopepa_network, runid);
    X_mRNA = importdata(loadfile);
    X = X_mRNA.data;
    
    % For the Biopepa data move around the light variable so it does not
    % get into our way of extracting the proper predictor data and response degradation term
    if(BIOPEPA_DATA) 
        % Put light*P (at column 3) as last entry, to make it easier for the parent-set
        % generation process
        X = [X, X(:,3)]; % Append the third feature, i.e. light*P for Biopepa to the end
        light_var = X(:,3);        % .. and save to additional use as indicator
        X(:,3) = [];               % delete the 3rd feature, since we moved it to the last column
    end
    
    
    % The feature response node is the degradation term, get it here:
    X_degrad = X(:,response_node)';
    
    % load the response (gradient) data
    loadfile = sprintf(['Data/biopepa_%s_%s_id%i.csv'], biopepa_network, gradient, runid);
    Y = importdata(loadfile);
    Y = Y.data;
    y = Y(:,response_node);
 
    
    %
    % Biopepa specific preprocessing of the data. We have prior knowledge
    % that specific observations can be excluded because the corresponding 
    % response node was knocked-down for the experiment.
    %
    % E.g. if response node 1 was selected, we know that in the samples 1:13
    %   the mRNA for node 1 was knocked-down. As an effect the gradient for
    %   the response node is flat and close to zero. Thus it makes sense to
    %   exclude these samples. 
    %
    % We have checked that the inference improves if we take this step.
    %
    if(BIOPEPA_DATA) 
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

        % This is Biopepa specific stuff. 
        % Take out observations if the current response is one of the below and if this feature is enabled
        if(response_node == 1) exclude_obs = 1:13;  end
        if(response_node == 2) exclude_obs = 14:26; end
        if(response_node == 6) exclude_obs = 27:39; end
        if(response_node == 4) exclude_obs = 40:52; end
        if(response_node == 5) exclude_obs = 40:52; end
        
        nr_obs = size(X,1);

        % Only keep samples which are not in the exclusion list
        take_these_obs = setdiff(1:nr_obs, exclude_obs);

        % Extract those samples defined above 
        X_degrad = X_degrad(take_these_obs);
        y = y(take_these_obs);
        X = X(take_these_obs,:);
    
    end
    
    
    
        
    % Pack it together into one data structure.
    Data.X = X';
    Data.X_degrad = X_degrad;
    Data.y = y;
    if(BIOPEPA_DATA)
        Data.light_var = light_var;
    end

end

