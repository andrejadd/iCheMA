
# iCheMA

This repository contains the Matlab implementation of iCheMA with adaptation to the Biopepa data set, which includes mRNA and gradient measurements from a stochastic process simulating transcriptional regulation. iCheMA was described in detail here: 

Aderhold, A., Grzegorczyk, M., and Husmeier, D. (2016). Approximate Bayesian inference in semi-mechanistic models. Statistics and Computing, 1-38.

iCheMA is a modified variant of CheMA introduced in 

Oates, C.J., Dondelinger, F., Bayani, N., Korkola, J., Gray, J.W.,
Mukherjee, S.: Causal network inference using biochemical kinetics.
Bioinformatics 30(17), i468¡Vi474 (2014)


## Overview

This folder contains the following content:

     * main_ICHEMA.m with function main_ICHEMA(DATA, node, parents, inhibition_vec, n_iterations)
       This function takes the following arguments:
       	    ** DATA - contains the design matrix with all features (DATA.X), the response gradient (DATA.y) and the degradation terms (DATA.X_degrad)
	    ** node - the response node for which to calculate the log likelihood
	    ** parents - a vector of indices that identify the predictor variables in DATA.X, can be empty, i.e. [] with no predictor
	    ** inhibition_vec - a matrix with 0 (activation) and 1 (inhibition) flags that define the type of term of each predictor.
	       each row corresponds to one setup
	    ** n_iterations - number of MCMC iterations.

     * run_Example_Biopepa.m 
       This function reads in the Biopepa data, generates the parent sets and the inhibition matrix that is needed for main_ICHEMA().
       Use this script as a template how to run the method. Modify load_Data() and getParentSets() to you own needs. The script will
       save the marginal log likelihood for each response and corresponding parent set configurations including the different activation
       and inhibition setups into the directory 'Results/[Gradient]' where [Gradient] specifies the type of Gradient. For Biopepa it can be
       the RBFGradient (analytic) or coarseGradient (numerical). This is specified inside run_Example_Biopepa.m

     * evaluate_Results.m.
       Reads out all the results from the directory 'Results/[Gradient]' and calculates the posterior probability for each response
       and predictor pair. In addition it will calculate an AUROC and AUPREC score in the case that a gold standard is provided. For Biopepa,
       the gold standard definitions are stored in the directory 'Data/Goldstandards'. The final results are saved to 'Results/EVAL'.

     * Directory 'Scripts' includes all the remaining scripts required by main_ICHEMA() and evaluate_Results().

     * Directory 'Data' contains the concentration and gradient data. 


## Usage

1. execute run_Example_Biopepa.m

This script runs iCheMA on several network and instances of the Biopepa data. The Biopepa data requires a couple of special preprocessing steps, which are implemented in load_Data(). This function takes a flag 'BIOPEPA_DATA', which can be set to zero in order to turn Biopepa specific processing off. Another important issue is the generation of different parent set configurations. The function getParentSets() builds the parent and activator/inhibitor combinations given a maximal fan-in. Furthermore, the number of available predictor variables influences the amount of parent sets, which has to be taken into account when modifying the code for new data (see getParentSets()).

Note: Running iCheMA on a single data set can take some time depending on the number of MCMC iterations, number of predictors, etc.. 

2.  execute evalute_Results.m

The script will read in the previously calculated marginal log likelihoods of all the parent set configurations and corresponding activation and inhibition term setups. It will produce a posterior probability score for each response and predictor variable given these scores. Finally, if the gold standard is known, as it is the case for the Biopepa data, the area under the ROC curve (AUROC) and the area under the PREC curve are calculated. These values are stored in the directory 'Results/EVAL'. 

Modify this script according to your needs. 
     
