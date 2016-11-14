

function [SAMPLE] = MCMC_SAMPLE_REDUCED_FIX_K(node,parents,inhibition_vec,prod_ind,DATA,n_iterations,K_star)

    a_SNR = 2;
    b_SNR = 0.2;

    y = DATA.y;

    nu_sig2 = 0.01;

    if(prod_ind==0)
        n_parents = length(parents);
    else
        n_parents = 1; 
    end

    n_obs  = length(y);
    V      = ones(n_parents+1,1); 
    sigma  = 1;
    delta  = 1;
    sigma2 = sigma^2;

    % Some fixed hyperparameters:
    nu = 0.5;

    % Store the initial values:
    SAMPLE.K{1}       = K_star;
    SAMPLE.V{1}       = V;
    SAMPLE.sigma{1}   = sigma;
    SAMPLE.delta{1}   = delta;
    SAMPLE.log_LL     = zeros(n_iterations,1);
    SAMPLE.log_Score  = zeros(n_iterations,1);


    % Build the design matrix:
    [D] = DESIGN_MATRIX_BUILDER(DATA, parents, inhibition_vec, prod_ind, K_star);
    % and compute the residuals:
    res = y - D*V;

    log_LL    = sum(lognormpdf_new(res,0,sigma));
    log_Score = log_LL + sum(lognormpdf_new(K_star,1,sqrt(nu))) + sum(lognormpdf_new(V,1,sqrt(delta*sigma^2))) + log(gampdf(sigma2^(-1),nu_sig2/2,(nu_sig2/2)^(-1)))  - sum(log(1-normcdf(zeros(length(K_star),1),1,sqrt(nu)))) - sum(log(1-normcdf(zeros(n_parents+1,1),1,sqrt(delta*sigma^2))));    

    SAMPLE.log_LL(1)    = log_LL;
    SAMPLE.log_Score(1) = log_Score;

    for t=2:(n_iterations+1)

        % Step1: 
        % Sample sigma and V via a Gibbs Sampling step

        [D] = DESIGN_MATRIX_BUILDER(DATA,parents,inhibition_vec,prod_ind,K_star);

        % The design matrix D is a (n_obs)-by-(n_parents+1) matrix 

        % Sufficient statistics for the regression parameters:
        % P(V|delta,sigma^2,(y,K,tau)) = N(mu_gibbs,Sigma_gibbs)

        SIGMA = inv(delta^(-1) * eye(n_parents+1) + D'*D);

        mu_gibbs    = SIGMA * (delta^(-1) * ones(n_parents+1,1) + D'*y);
        Sigma_gibbs = sigma^2 * SIGMA;

        % Re-sample the regression parameter vector V:

        OKAY    = 0;
        Counter = 0;

        while((OKAY==0) && (Counter <= 10))

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
        % P(sigma^2|delta,y,K) = IG(a,b)

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


        % No Step2, as K_star is kept fixed 

        [D_star] = DESIGN_MATRIX_BUILDER(DATA, parents,inhibition_vec,prod_ind,K_star);
        res_star = y - D_star*V;

        log_LL_star    = sum(log(normpdf(res_star,0,sigma))); 
        log_Score_star = log_LL_star + sum(lognormpdf_new(K_star,1,sqrt(nu))) + sum(lognormpdf_new(V,1,sqrt(delta*sigma^2))) + log(gampdf(sigma2^(-1),nu_sig2/2,(nu_sig2/2)^(-1)))  - sum(log(1-normcdf(zeros(length(K_star),1),1,sqrt(nu)))) - sum(log(1-normcdf(zeros(n_parents+1,1),1,sqrt(delta*sigma^2)))); 

        SAMPLE.log_LL(t)    = log_LL_star;
        SAMPLE.log_Score(t) = log_Score_star;  
        SAMPLE.K{t}         = K_star;

    end

return

