
function [parent_set, inhibitors] = getParentSets(parent_nodes, parent_set_id, max_fanin)
    

 % This counter is used to identify the parent set generated with
    % recombination of different parent ids for different parent sizes
   
    parent_set = [];  % the empty set corresponds to parent_set_id = 1

    if parent_set_id > 1

        % start with 2, since p_counter == 1 is empty set []
        p_counter = 2;

        found_parent_set = 0;
        
        for ii = 1:max_fanin          

            comb_tmp = nchoosek(parent_nodes, ii);                                                                                                                     

            for j = 1:size(comb_tmp,1) 

                % if the 
                if p_counter == parent_set_id;
             
                    parent_set = comb_tmp(j,:);
                    found_parent_set = 1;
                    break;
                end
                
                % increase parent counter
                p_counter = p_counter + 1;
            end
            
            if found_parent_set
                break;
            end
            
        end
    end                                                                                                                                                                     
    
    % check if parent_set was selected inside boundary
    if (parent_set_id > 1) && isempty(parent_set)
         err = MException('user:badvalue', ' parent_set_id index out of bounds.');
         throw(err);
    end
    
    fprintf('selected parent set for id %i:', parent_set_id);
    parent_set

    
    
    nr_parents = length(parent_set);

    fprintf('number of parents (including the degradation term): %i\n', (nr_parents+1));
    
    % create all possible combinations of activators/inhibitors as a binary
    % vector where 0 --> activator, 1 --> inhibitor 
    inhibitors = [];
    
    % if the parent set is empty, do not try to recombine binaries
    if ~isempty(parent_set)
        fact_vec = ones(1,nr_parents) + 1;
        inhibitors = fullfact(fact_vec) - 1;
    end
        
    
    
end