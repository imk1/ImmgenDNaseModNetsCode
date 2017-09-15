function obj = calculate_objectiveDNase ( modmean, tf_expression, assign, nonEmptyModules, W, DNaseRel, cellTypeMap, DNasePairWeights, e )

% Code written by Su-In Lee
% Modified by Irene Kaplow

% Calculates the full objective for DNase Lirnet

% input:
%       modmean: the average expression of the genes in each module for
%           each cell type -- matrix is num_mods x n, where num_mods is the
%           number of modules and n is the number of cell types
%       tf_expression:  tf expression, size rxn where r is the number of
%           regulators, each TF's expression has been normalized (if necessary)
%       assign: module assignments
%       nonEmptyModules: list of modules that contain at least 1 gene
%       DNaseRel: the relevance of each regulator for each gene in each
%           cell type -- matrix that is r x g x i, where r is the number of
%           regulators, g is the number of genes, and i is the number of cell
%           types with DNase
%       cellTypeMap: maps each cell type to the closest cell type with
%           near-by DNase -- vector of length n, where n is the number of 
%           cell types and where each entry has the index of the closest
%           cell type with DNase
%       DNasePairWeights: initial values of meta-feature weights
%       num_mod: is the number of modules
%       e : penalty for l2 regularization for DNase feature weights
% output:
%       obj: the value of the objective

DNaseRelAverage = zeros(size(W',1), size(W',2), size(DNaseRel,3));
for k = 1:length(nonEmptyModules)
    % Compute the average relevance for each regulator across genes in the
    % module
    geneLocationsIndexes = find(ismember(assign, nonEmptyModules(k)));
    for i = 1:size(DNaseRel, 3)
        % Iterate through cell types and find the average regulator
        % relevance for each
        DNaseRelAverage(:,nonEmptyModules(k),i) = mean(DNaseRel(:,geneLocationsIndexes,i), 2)'; % Using 1 l1 penalty that is the average of all of the l1 penalties for the genes in the module
    end
end

weightTimesDNaseRelAverage = zeros(length(nonEmptyModules), size(DNaseRelAverage, 1), size(DNaseRelAverage, 3));
for i = 1:size(DNaseRel, 3)
    % Iterate through cell types and multiply the absolute value of the 
    % weight by the meta-prior for each cell type
    weightTimesDNaseRelAverage(:,:,i) = abs(W(nonEmptyModules,:)) .* (DNaseRelAverage(:,nonEmptyModules, i)');
end

mainObj = .5*sum( (modmean(nonEmptyModules,:) - W(nonEmptyModules,:) * tf_expression).^2, 2);

objectiveWithMetapriorVarsCellType = zeros(length(nonEmptyModules), size(DNaseRelAverage,1));
for ct = 1:length(cellTypeMap)
    % Iterate through all cell types and add the information for the prior
    % for the nearest cell type with DNase to the inner part of the
    % objective
    objectiveWithMetapriorVarsCellType = objectiveWithMetapriorVarsCellType + weightTimesDNaseRelAverage(:,:,cellTypeMap(ct)) - permute(log(DNaseRelAverage(:,nonEmptyModules,cellTypeMap(ct))), [2 1 3]);
end

obj = sum(mainObj + sum(objectiveWithMetapriorVarsCellType, 2)) + (e*sum(DNasePairWeights.^2));