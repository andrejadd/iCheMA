
function [log_LL_marginal] = main_ICHEMA(DATA, node, parents, inhibition_vec, n_iterations)

 
    % Initialize some parameters
    factor = 0.1;
    nu_sig2 = 0.01;
    a_SNR   = 2;
    b_SNR   = 0.2;
    nu      = 0.5;
    prod_ind = 0;
    y = DATA.y;
    n_obs = length(y);

    
    % Perform a full MCMC simulation
    [SAMPLE] = MCMC_SAMPLE(node,parents,inhibition_vec,prod_ind, DATA, n_iterations, factor);

    % Use the MCMC SAMPLING PHASE trajectory to determine a high-density point
    % after a burn_in of half the iteration length:
    n_burn = floor(n_iterations/2);
    [~,index] = max(SAMPLE.log_Score(n_burn + 1:end));
    index       = index + n_burn;

    K_star      = SAMPLE.K{index};
    V_star      = SAMPLE.V{index};
    sigma_star  = SAMPLE.sigma{index};
    delta_star  = SAMPLE.delta{index};


    % extract samples and calculate autocorrelation
    M = n_iterations - n_burn;

    for jj = 1:length(SAMPLE.V{1})
        ss1 = [];

        for ii = n_burn:n_iterations
            ss1 = [ss1, SAMPLE.V{ii}(jj)];
        end
        
        sample_autoc = autocorr(ss1);
        
        ESS = M / (1 + 2 * sum(sample_autoc));

        % fprintf('sample %i, ESS=%d (sum p = %d)\n', jj, ESS, sum(sample_autoc));  


    end



    % Now perform a reduced MCMC simulations WHERE 
    % K_star is kept fixed...
    [SAMPLE2] = MCMC_SAMPLE_REDUCED_FIX_K(node, parents, inhibition_vec, prod_ind, DATA, n_iterations, K_star);

    %%%%%%%%%%% NOT REQUIRED
    % Now perform a reduced MCMC simulation WHERE K_star and sigma_star are kept fixed...
    [SAMPLE3] = MCMC_SAMPLE_REDUCED_FIX_K_sigma(node, parents, inhibition_vec, prod_ind, DATA, n_iterations, K_star, sigma_star);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if(prod_ind==0)
        n_parents = length(parents);
    else
        n_parents = 1; 
    end

  

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Generate a 'variate' for SAMPLE2, where K=K_star was kept fixed

    for j=1:n_iterations
        x            = factor*randn(length(parents),1);
        SAMPLE2.K{j} = abs(K_star + x);   
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Compute the (reduced) 'marginal ordinates', i.e. compute Eq (17)
    % \pi(K_star|delta,sigma,V,tau)

    numerator   = 0;
    denominator = 0;

    n_burn = floor(n_iterations/2);

    % NUMERATOR:
    for g=(n_burn+1):n_iterations

        % LL_star and LL_g
        [D_star] = DESIGN_MATRIX_BUILDER(DATA, parents,inhibition_vec,prod_ind,K_star    );
        [D_g]    = DESIGN_MATRIX_BUILDER(DATA, parents,inhibition_vec,prod_ind,SAMPLE.K{g});

        res_star = y - D_star*SAMPLE.V{g};
        res_g    = y - D_g   *SAMPLE.V{g};

        log_LL_star = sum(lognormpdf_new(res_star,0,SAMPLE.sigma{g})) + sum(lognormpdf_new(K_star,     1,sqrt(nu)))  - sum(log(1-normcdf(zeros(length(K_star),1),1,sqrt(nu))));
        log_LL_g    = sum(lognormpdf_new(res_g,   0,SAMPLE.sigma{g})) + sum(lognormpdf_new(SAMPLE.K{g},1,sqrt(nu)))  - sum(log(1-normcdf(zeros(length(K_star),1),1,sqrt(nu)))) ;

        A_g_star = min(exp(log_LL_star - log_LL_g),1);

        numerator = numerator + prod(normpdf(K_star - SAMPLE.K{g},0,factor)) * A_g_star; % factor is the standard-deviation
    end

    % DENOMINATOR:
    for j=(n_burn+1):n_iterations

        % LL_star and LL_g
        [D_star] = DESIGN_MATRIX_BUILDER(DATA,parents,inhibition_vec,prod_ind,K_star     );
        [D_j]    = DESIGN_MATRIX_BUILDER(DATA,parents,inhibition_vec,prod_ind,SAMPLE2.K{j});

        res_star = y - D_star*SAMPLE2.V{j};
        res_j    = y - D_j   *SAMPLE2.V{j};

        log_LL_star = sum(lognormpdf_new(res_star,0,SAMPLE2.sigma{j})) + sum(lognormpdf_new(K_star,     1,sqrt(nu)))  - sum(log(1-normcdf(zeros(length(K_star)      ,1),1,sqrt(nu))));
        log_LL_j    = sum(lognormpdf_new(res_j,   0,SAMPLE2.sigma{j})) + sum(lognormpdf_new(SAMPLE2.K{j},1,sqrt(nu))) - sum(log(1-normcdf(zeros(length(SAMPLE2.K{j}),1),1,sqrt(nu))));

        A_star_j = min(exp(log_LL_j - log_LL_star),1);

        denominator = denominator + A_star_j;
    end

    log_marginal_ordinate_1 = log(numerator) - log(denominator);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Compute the three remaing (reduced) 'marginal ordinates', i.e. compute Eq (17)
    % Note that sigma, delta and V are sampled by Gibbs-sampling steps

    % P(sigma_star|K_star)
    numerator = 0;
    n_burn    = floor(n_iterations/2);

    % ONLY NUMERATOR:
    for g=(n_burn+1):n_iterations

        [D] = DESIGN_MATRIX_BUILDER(DATA,parents,inhibition_vec,prod_ind,K_star);

        % Sufficient statistics for sigma^2, where V has been integrated out
        % P(sigma^2|delta,K_star,tau) = IG(a,b)

        a = nu_sig2/2 + 0.5 * (n_obs+n_parents+1);
        b = nu_sig2/2 + 0.5 * ((y-D*SAMPLE2.V{g})' * (y-D*SAMPLE2.V{g}) + SAMPLE2.delta{g}^(-1) * (SAMPLE2.V{g}-ones(n_parents+1,1))' *  (SAMPLE2.V{g}-ones(n_parents+1,1)));

        numerator = numerator + gampdf(sigma_star^(-2),a,1/b);

    end

    log_marginal_ordinate_2 = log(numerator) - log(n_iterations-n_burn);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % P(delta_star|K_star,sigma_star)
    numerator   = 0;
    n_burn = floor(n_iterations/2);

    % ONLY NUMERATOR:
    for g=(n_burn+1):n_iterations

        [D] = DESIGN_MATRIX_BUILDER(DATA,parents,inhibition_vec,prod_ind,K_star);

        % Sufficient statistics for the regression parameters:
        % P(V|delta,sigma,y,K) = N(mu_gibbs,Sigma_gibbs)

        a = a_SNR + 0.5 * (n_parents+1);
        b = b_SNR + 0.5 * sigma_star^(-2) * (SAMPLE3.V{g}-ones(n_parents+1,1))' * (SAMPLE3.V{g}-ones(n_parents+1,1));

        numerator = numerator + gampdf(delta_star^(-1),a,1/b);

    end

    log_marginal_ordinate_3 = log(numerator) - log(n_iterations-n_burn);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % P(V_star|K_star,sigma_star,delta_star)

    [D] = DESIGN_MATRIX_BUILDER(DATA,parents,inhibition_vec,prod_ind,K_star);


    SIGMA =  inv(delta_star^(-1) *eye(n_parents+1) + D'*D);

    mu_gibbs    = SIGMA * (delta_star^(-1) * ones(n_parents+1,1) + D'*y);
    Sigma_gibbs = sigma_star^2 * SIGMA;

    log_marginal_ordinate_4 = log(mvnpdf(V_star,mu_gibbs,Sigma_gibbs)) - log(mvncdf(zeros(n_parents+1,1),Inf*ones(n_parents+1,1),mu_gibbs,Sigma_gibbs));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Compute the log of marginal likelihood and prior 
    % for the parameters (K_star,V_star,sigma_star,delta_star,tau_star)

    [D_star] = DESIGN_MATRIX_BUILDER(DATA,parents,inhibition_vec,prod_ind,K_star);
    res_star = y - D_star*V_star;

    log_LL           = sum(lognormpdf_new(res_star,0,sigma_star)); 
    log_LL_and_prior = log_LL + sum(lognormpdf_new(K_star,1,sqrt(nu))) - sum(log(1-normcdf(zeros(length(K_star),1),1,sqrt(nu)))) + log(gampdf(sigma_star^(-2),nu_sig2/2,(nu_sig2/2)^(-1)));
    log_LL_and_prior = log_LL_and_prior  + sum(lognormpdf_new(V_star,1,sqrt(delta_star*sigma_star^2)))- sum(log(1-normcdf(zeros(length(V_star),1),1,sqrt(delta_star*sigma_star^2))));

    %%%%%%%%%%%%%

    log_LL_marginal = log_LL_and_prior - (log_marginal_ordinate_1 + log_marginal_ordinate_2 + log_marginal_ordinate_3 + log_marginal_ordinate_4);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



