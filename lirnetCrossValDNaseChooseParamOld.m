function crossValParams = lirnetCrossValDNaseChooseParam(expression, tf_expression, assign, DNasePairFeats, cellTypeMap, DNasePairWeights, num_mod, maxiters, reclust, l1UpperList, l1LowerList, e, crossValFold)

% Run Lirnet with k-fold cross-validation and DNase-related priors that are 
% gene-specific

% input:
%       expression: the expression file, size gxn where g is the number of genes
%               and n is the number of cell types, each gene's expression has
%               been normalized (if necessary)
%       tf_expression:  tf expression, size rxn where r is the number of
%           regulators, each TF's expression has been normalized (if necessary)
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
%       DNasePairWeights: initial values of meta-feature weights
%       num_mod: is the number of modules
%       maxiters: the maximum number of iterations for the metaprior loop
%       reclust: number of times to move genes to better modules and repeat
%           the Lirnet process
%       l1UpperList: an array with the list of the upper bounds on the l1
%           penalty that will be tried
%       l1LowerList: an array with the list of the lower bounds on the l1
%           penalty that will be tried
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
%                   types
%               DNasePairWeights: matrix of size hp, where hp is the number 
%                   of DNase pairwise features, that contains the learned 
%                   metapriors for the DNase pairwise features
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
%           l1Upper: C0 value that gave the lowest validation error
%           l1Lower: C1 value that gave the lowest validation error


crossValParamsListPlus = {};
meanWeightedMeanValErrorsPlus = zeros(length(l1LowerList), 1);

for j = 1:length(l1LowerList)
    % Iterate through C1 values and try Lirnet with each
    fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
    fprintf('   C1 = %g, %g/ %g\n', l1LowerList(j), j, length(l1LowerList));
    fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n')
    crossValParamsList = {};
    meanWeightedMeanValErrors = [];
    l1UpperUsed = [];
    l1Lower = l1LowerList(j);
    l1UpperCount = 0;
    for i = 1:length(l1UpperList)
        % Iterate through C0 values and try Lirnet with each
        if l1Lower > l1UpperList(i)
            % C1 should not be > C0, so do not try this C0
            continue
        end
        l1UpperCount = l1UpperCount + 1;
        fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
        fprintf('   C0 = %g, %g/ %g\n', l1UpperList(i), i, length(l1UpperList));
        fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n')
        l1UpperUsed(l1UpperCount) = l1UpperList(i);
        l1Bounds = [l1UpperList(i), l1Lower];
        crossValParamsList{l1UpperCount} = lirnetCrossValDNase(expression, tf_expression, assign, DNasePairFeats, cellTypeMap, DNasePairWeights, num_mod, maxiters, reclust, l1Bounds, e, crossValFold);
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