function crossValParams = lirnetCrossValDNase(expression, tf_expression, assign, DNasePairFeats, cellTypeMap, DNasePairWeights, num_mod, maxiters, reclust, l1Bounds, e, crossValFold)

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
%       l1Bounds: an array with 2 entries: the upper bound on the l1
%           penalty followed by the lower bound on the l1 penalty
%       e : penalty for l2 regularization for DNase feature weights
%       crossValFold: cross-validation fold
% output:
%       crossValParams is a structure with the following fields:
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

numCellTypes = size(expression, 2);
valSetSize = round(numCellTypes/crossValFold);
randList = randperm(numCellTypes);

crossValParams.paramsList = {};
crossValParams.moduleSizes = zeros(crossValFold, num_mod);
crossValParams.numModules = zeros(crossValFold, 1);
crossValParams.numNonZeroWeights = zeros(crossValFold, num_mod);
crossValParams.meanNumNonZeroWeights = zeros(crossValFold, 1);
crossValParams.valErrors = zeros(crossValFold, num_mod);
crossValParams.meanValErrors = zeros(crossValFold, 1);

for i = 1:crossValFold
    % Run Lirnet with cross-validation crossValFold times
    
    fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
    fprintf('   CROSSVALITER = %g/ %g\n', i, crossValFold);
    fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
    
    valSetIndexes = randList((i-1)*valSetSize+1:min([numCellTypes; i*valSetSize]));
    trainSetIndexes = setdiff(randList, valSetIndexes);
    expressionTrain = expression(:,trainSetIndexes);
    tf_expressionTrain = tf_expression(:,trainSetIndexes);
    
    params = lirnetDNase(expressionTrain, tf_expressionTrain, assign, DNasePairFeats, cellTypeMap, DNasePairWeights, num_mod, maxiters, reclust, l1Bounds, e);
    crossValParams.paramsList{i} = params;
    
    moduleSizes = zeros(num_mod, 1);
    for j = 1:num_mod
        % Iterate through modules and find the size of each module
        for k = 1:size(params.assign, 1)
            % Iterate through genes and find the number in each module
            if params.assign(k, size(params.assign, 2)) == j
                % The current gene is in the current module, so add 1 to
                % the size of the current module
                moduleSizes(j) = moduleSizes(j) + 1;
            end
        end
    end
    crossValParams.moduleSizes(i,:) = moduleSizes';
    crossValParams.numModules(i) = length(find(moduleSizes > 0));
    
    numNonZeroWeights = zeros(num_mod, 1);
    for j = 1:num_mod
        % Iterate through modules and find the number of TFs in each module
        % with non-0 weights
        numNonZeroWeights(j) = length(find(params.weights{length(params.weights)}(j,:) ~= 0));
    end
    crossValParams.numNonZeroWeights(i,:) = numNonZeroWeights';
    crossValParams.meanNumNonZeroWeights(i) = mean(numNonZeroWeights(find(moduleSizes > 0)));
    
    expressionVal = expression(:, valSetIndexes);
    tf_expressionVal = tf_expression(:, valSetIndexes);
    
    valErrors = zeros(num_mod, 1);
    for j = 1:num_mod
        % Iterate through modules and compute the regression error for each
        % module using the validation set
        g = find(params.assign(:,size(params.assign, 2))==j);
        if isempty(g)
            % The module is empty, so do not consider it
            continue
        end
        pred = tf_expressionVal'*params.weights{length(params.weights)}(j,:)'; % predictions
        valErrors(j) = sum ( sum ( (expressionVal(g,:) -  repmat(pred',length(g),1)).^2)); %sse    
    end
    crossValParams.valErrors(i,:) = valErrors';
    crossValParams.meanValErrors(i) = mean(valErrors(find(moduleSizes > 0)));
    crossValParams.weightedMeanValErrors(i) = (valErrors' * moduleSizes)/sum(moduleSizes);
end
