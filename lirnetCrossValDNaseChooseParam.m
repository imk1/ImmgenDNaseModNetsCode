function crossValParams = lirnetCrossValDNaseChooseParam(expression, tf_expression, standReq, assign, DNasePairFeats, cellTypeMap, DNaseRel, num_mod, maxiters, reclust, l1UpperList, l1LowerList, d, e, crossValFold)

% Run Lirnet with k-fold cross-validation and DNase-related priors that are 
% gene-specific

% input:
%       expression: the expression file, size gxn where g is the number of genes
%               and n is the number of cell types, each gene's expression has
%               been normalized (if necessary)
%       tf_expression:  tf expression, size rxn where r is the number of
%           regulators, each TF's expression has been normalized (if necessary)
%       standReq: indicates what standardization has been done -- -1 means
%           that no standardization has been done, 0 means that the mean
%           has been subtracted but division by the standard deviation has
%           not been done, and 1 means that full standardization has been
%           done
%       assign: initial module assignments, leave it empty for k-means
%       clustering
%       DNasePairFeats: pairwise metafeatures for DNase sites -- cell array
%           where each cell has a different pairwise feature for a gene and 
%           a TF, each cell's entry is g x r x i, where g is the number of
%           genes, r is the number of regulators (assumes that all 
%           features relating to a DNase site that are not near a gene have 
%           value 0 and that PWM scores have already been incorporated into
%           features), and i is the number of cell types with DNase
%       cellTypeMap: maps each cell type to the closest cell type with
%           near-by DNase -- vector of length n, where n is the number of 
%           cell types and where each entry has the index of the closest
%           cell type with DNase
%       DNaseRel: the initial values of relevance of each regulator for 
%           each gene in each cell type -- matrix that is r x g x n, where 
%           r is the number of regulators, g is the number of genes, and n 
%           is the number of cell types with DNase
%       num_mod: is the number of modules
%       maxiters: the maximum number of iterations for the metaprior loop
%       reclust: number of times to move genes to better modules and repeat
%           the Lirnet process
%       l1UpperList: an array with the list of the upper bounds on the l1
%           penalty that will be tried
%       l1LowerList: an array with the list of the lower bounds on the l1
%           penalty that will be tried
%       d: L2 parameter for feature weights
%       e : penalty for l2 regularization for DNase feature weights
%       crossValFold: cross-validation fold
% output:
%       crossValParams is a structure with the following fields and is 
%           and is chosen to be the one for the C0 and C1 pair with the 
%           lowest average validation error:
%           paramsList: a cell array in which each entry has the parameters
%               from a round of cross-validation, where the parameters are 
%               in a params structure, which has several relevant fields:
%               W0: matrix of size num_mods x r, where are is the number of
%                   regulators, with the initial learned weights for each
%                   regulator
%               errors: is the sum squared error
%               assign: has reclust columns, each column is the assigment 
%                   of genes into cluster for the given iteration (reclust 
%                   iterations in total) (note: do not use more than 2 
%                   iterations if the penalty for l1 regularization is 
%                   stringent)
%               weights: cell array in which each entry is a matrix of size 
%                   num_mods x r, coefficients for each module, where the 
%                   number of entires is equal to the number of 
%                   re-clustering events
%               DNaseRel: the relevance of each regulator for each gene in each
%                   cell type -- matrix that is r x g x n, where r is the number of
%                   regulators, g is the number of genes, and n is the number of cell
%                   types with DNase
%               DNasePairWeights: cell array in which each entry is a matrix of size hp x r, where hp is the number of
%                   DNase pairwise features and r is the number of regulators, that contains the learned
%                   metapriors for the DNase pairwise features, and there is an
%                   entry for each re-clustering iteration
%           moduleSizes: a crossValFold x num_mod array, where each entry
%               (i,j) is the size of final module j in cross-validation
%               iteration i
%           numModules: an array of size crossValFold, where each entry is
%               the number of non-empty modules in the corresponding
%               cross-validation iteration
%           numNonZeroWeights: a crossValFold x num_mod array, where each
%               entry (i,j) is the number of TFs with non-zero weights for
%               module j in cross-validation iteration i
%           meanNumNonZeroWeights: an array of size crossValFold, where
%               each entry is the average of the numbers of TFs with
%               non-zero weights across all non-empty modules in the
%               corresponding iteration
%           valErrors: a crossValFold x num_mod array, where each entry
%               (i,j) is the sse on the validation set of genes in module j
%               in cross-validation iteration i
%           meanValErrors: an array of size crossValFold, where each entry
%               is the average of the sses across all non-empty modules in
%               the corresponding iteration
%           weightedMeanValErrors: an array of size crossValFold, where
%               each entry is the weighted average of the sses across the
%               modules in the corresponding iteration, where weightings
%               are done based on module size
%           l1Upper: L1 penalty upper boound value that gave the lowest 
%               validation error
%           l1Lower: L1 penalty lower bound that gave the lowest validation 
%   `           error
%           valCorrs: correlation between predicted expression and true
%               expression for each gene for each validation set, array is
%               size crossValFold x g, where g is the number of genes
%           meanValCorrs: array of size crossValFold, where each entry
%               is the average correlation across genes of gene expression
%               and predicted gene expression for that cross-validation
%               iteration


crossValParamsListPlus = {};
meanWeightedMeanValErrorsPlus = zeros(length(l1LowerList), 1);

for j = 1:length(l1LowerList)
    % Iterate through L1 penalty lower bound values values and try Lirnet with each
    fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
    fprintf('   C1 = %g, %g/ %g\n', l1LowerList(j), j, length(l1LowerList));
    fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n')
    crossValParamsList = {};
    meanWeightedMeanValErrors = [];
    l1UpperUsed = [];
    l1Lower = l1LowerList(j);
    l1UpperCount = 0;
    for i = 1:length(l1UpperList)
        % Iterate through L1 penalty upper bound values and try Lirnet with each
        l1UpperCount = l1UpperCount + 1;
        fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
        fprintf('   C0 = %g, %g/ %g\n', l1UpperList(i), i, length(l1UpperList));
        fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n')
        l1UpperUsed(l1UpperCount) = l1UpperList(i);
        l1Bounds = [l1UpperList(i), l1Lower];
        crossValParamsList{l1UpperCount} = lirnetCrossValDNase(expression, tf_expression, standReq, assign, DNasePairFeats, cellTypeMap, DNaseRel, num_mod, maxiters, reclust, l1Bounds, d, e, crossValFold);
        meanWeightedMeanValErrors(l1UpperCount) = mean(crossValParamsList{l1UpperCount}.weightedMeanValErrors);
    end
    [~, minIndex] = min(meanWeightedMeanValErrors);
    crossValParamsListPlus{j} = crossValParamsList{minIndex};
    crossValParamsListPlus{j}.l1Upper = l1UpperUsed(minIndex);
    meanWeightedMeanValErrorsPlus(j) = mean(crossValParamsListPlus{j}.weightedMeanValErrors);
end

[~, minIndexPlus] = min(meanWeightedMeanValErrorsPlus);
crossValParams = crossValParamsListPlus{minIndexPlus};
crossValParams.l1Lower = l1LowerList(minIndexPlus);