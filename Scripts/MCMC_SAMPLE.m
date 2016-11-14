

function [SAMPLE] = MCMC_SAMPLE(node, parents, inhibition_vec, prod_ind, DATA, n_iterations, factor)

    nu_sig2 = 0.01;
    a_SNR   = 2;
    b_SNR   = 0.2;
    
    % get response vector
    y = DATA.y;
    
    % No. of regulators (ignoring the degradation process)
    if(prod_ind == 0)
        n_parents = length(parents);
    else
        n_parents = 1;   
    end

    % No. of time points:
    n_obs = length(y);

    % Initialisation
    K      = ones(n_parents,1);
    V      = ones(n_parents+1,1); 
    sigma  = 1;       
    sigma2 = sigma^2; 
    delta  = 1;

    nu = 0.5; % hyperparameter

    % Store the initial values:
    SAMPLE.K{1}       = K;
    SAMPLE.V{1}       = V;
    SAMPLE.sigma{1}   = sigma;
    SAMPLE.delta{1}   = delta;
    SAMPLE.log_LL     = zeros(n_iterations,1);
    SAMPLE.log_Score  = zeros(n_iterations,1);

    % Build the design matrix:
    [D] = DESIGN_MATRIX_BUILDER(DATA, parents, inhibition_vec, prod_ind, K);

    % Compute the residuals:
    res = y - D*V;

    log_LL    = sum(lognormpdf_new(res,0,sigma)); % the log-likelihood
  
    % Add the log-prior to the log-likelihood 
    log_Score = log_LL + sum(lognormpdf_new(K,1,sqrt(nu))) - sum(log(1-normcdf(zeros(length(K),1),1,sqrt(nu)))) + sum(lognormpdf_new(V,1,sqrt(delta*sigma^2))) - sum(log(1-normcdf(zeros(n_parents+1,1),1,sqrt(delta*sigma^2)))) + log(gampdf(sigma2^(-1),nu_sig2/2,(nu_sig2/2)^(-1)));    

    SAMPLE.log_LL(1)    = log_LL;
    SAMPLE.log_Score(1) = log_Score;

    % Start the MCMC-simulation

    for t=2:(n_iterations+1)

        % Step1: 
        % Sample sigma and V via a Gibbs Sampling step

        [D] = DESIGN_MATRIX_BUILDER(DATA, parents, inhibition_vec, prod_ind, K);

        % Just note that the design matrix D is a (n_obs)-by-(n_parents+1) matrix 

        % Sufficient statistics for the regression parameters:
        % P(V|delta,sigma^2,(y,K,tau)) = N(mu_gibbs,Sigma_gibbs)

        SIGMA = inv(delta^(-1) * eye(n_parents+1) + D'*D);

        mu_gibbs    = SIGMA * (delta^(-1) * ones(n_parents+1,1) + D'*y);
        Sigma_gibbs = sigma^2 * SIGMA;

        % Re-sample the regression parameter vector V:

        % The multivariate Gaussian distribution of V is truncated to "V>=0"
        % So sample (up to 5-times) until this condition is met

        % As it appears that only V(1) tends to be negative:

        OKAY    = 0;
        Counter = 0;

        while(OKAY==0 && Counter <=10)

            V = mvnrnd(mu_gibbs,Sigma_gibbs)';

            if(min(V)>=0)
               OKAY=1;
            else
               Counter = Counter +1; 
            end    
        end


         if (Counter>=10)
             V = mvnrnd(mu_gibbs,Sigma_gibbs)';
             indicis = find(V<0);
             V(indicis) = 0;
         end


        % Sufficient statistics for sigma^2
        % P(sigma^2|delta,(y,K,tau)) = IG(a,b)

        a = nu_sig2/2 + 0.5 * (n_obs+length(V));
        b = nu_sig2/2 + 0.5 * (  (y-D*V)' * (y-D*V) + delta^(-1) * (V-ones(n_parents+1,1))' *  (V-ones(n_parents+1,1))  );

        % Re-sample sigma
        inv_sigma2 = gamrnd(a,1/b); % one has to invert b (???)
        sigma2     = inv_sigma2^(-1);
        sigma      = sqrt(sigma2);

        % Sufficient statistics for delta:
        % P(delta|V,sigma^2,(y,K,tau)) = IG(a,b)

        a = a_SNR + 0.5 * (n_parents+1);
        b = b_SNR + 0.5 * sigma^(-2) * (V-ones(n_parents+1,1))' * (V-ones(n_parents+1,1));

        % Re-sample delta: 
        inv_delta = gamrnd(a,1/b); 
        delta     = inv_delta^(-1);

        % Store the new sample of V, sigma and delta:
        SAMPLE.V{t}     = V;
        SAMPLE.sigma{t} = sigma;
        SAMPLE.delta{t} = delta;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Step2: 
        % Sample the vector K via a Metropolis Hastings step:

        % The current vector K:
        K_OLD = SAMPLE.K{t-1};
        x     = factor*randn(length(K_OLD),1);

        % Propose a new candidate for K:
        K_NEW = abs(K_OLD + x);

        % LL_old
        [D_old] = DESIGN_MATRIX_BUILDER(DATA, parents, inhibition_vec, prod_ind, K_OLD);
        res_old = y - D_old*V;

        % LL_new
        [D_new] = DESIGN_MATRIX_BUILDER(DATA, parents, inhibition_vec, prod_ind, K_NEW);
        res_new = y - D_new*V;

        % Compute the log-likelihoods...
        log_LL_old = sum(lognormpdf_new(res_old,0,sigma)); 
        log_LL_new = sum(lognormpdf_new(res_new,0,sigma));

        % ...and add the log-priors:
        log_Score_old = log_LL_old + sum(lognormpdf_new(K_OLD,1,sqrt(nu))) + sum(lognormpdf_new(V,1,sqrt(delta*sigma^2))) + log(gampdf(sigma2^(-1),nu_sig2/2,(nu_sig2/2)^(-1)))  - sum(log(1-normcdf(zeros(length(K_OLD),1),1,sqrt(nu)))) - sum(log(1-normcdf(zeros(n_parents+1,1),1,sqrt(delta*sigma^2)))); 
        log_Score_new = log_LL_new + sum(lognormpdf_new(K_NEW,1,sqrt(nu))) + sum(lognormpdf_new(V,1,sqrt(delta*sigma^2))) + log(gampdf(sigma2^(-1),nu_sig2/2,(nu_sig2/2)^(-1)))  - sum(log(1-normcdf(zeros(length(K_NEW),1),1,sqrt(nu)))) - sum(log(1-normcdf(zeros(n_parents+1,1),1,sqrt(delta*sigma^2)))); 

        log_inv_HR = 1; % the inverse Hastings Ratio is equal to one

        A = exp(log_Score_new - log_Score_old + log_inv_HR);

        u = rand(1);

        if (u<A) % accept K_new
            K = K_NEW;
            SAMPLE.log_LL(t)    = log_LL_new;
            SAMPLE.log_Score(t) = log_Score_new;
        else % or leave K unchanged, i.e. K=K_old
            K = K_OLD;
            SAMPLE.log_LL(t)    = log_LL_old;
            SAMPLE.log_Score(t) = log_Score_old;
        end

        % Store the new sample of K
        SAMPLE.K{t} = K;
    end

return






