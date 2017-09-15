function DNasePairWeights = learn_regul_metapriorDNase (DNasePairFeats, cellTypeMap, DNasePairWeights, W, assign, l1Upper, l1Lower, e)

% Code written by Su-In Lee
% Modified by Irene Kaplow

% Calculates the meta-priors for the DNase features

% input:
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
%       DNasePairWeights: values of meta-feature weights
%       W: a matrix of size num_mods x r, that contains the 
%           coefficients for each module, regulator pair
%       assign: module assignments
%       l1Upper: the upper bound on the l1 penalty (C0 in Lirnet paper)
%       l1Lower: the lower bound on the l1 penalty (C1 in Lirnet paper)
%       num_mod: is the number of modules
%       e : penalty for l2 regularization for DNase feature weights
% output:
%       DNasePairWeights: new values of meta-feature weights

fprintf('BOF::learn_regul_metaprior\n');

[DNasePairWeights,f,funEvals] = minConF_PQNDNase(@computeObjectiveWithMetapriorVars, DNasePairWeights, @projectOntoNegative, DNasePairFeats, cellTypeMap, W, assign, l1Upper, l1Lower, e);

fprintf('EOF::learn_marker_metaprior\n\n');
end


function [objectiveWithMetapriorVars, pgrad] = computeObjectiveWithMetapriorVars(DNasePairWeights, DNasePairFeats, cellTypeMap, W, assign, l1Upper, l1Lower, e)
% Computes the objective given a set of feature weights
u = abs(W)';
ngenes = size( DNasePairFeats{1},1 );
nCellTypes = size(DNasePairFeats{1}, 3);

DNaseRel = calculate_regul_relevanceDNase (l1Upper, l1Lower, DNasePairFeats, DNasePairWeights);
DNaseRelAverage = zeros(size(u,1), size(u,2), nCellTypes);
nonEmptyModules = [];
for k = 1:size(u, 2)
    % Compute the average relevance for each regulator across genes in the
    % module
    geneLocationsIndexes = find(ismember(assign, k));
    if isempty(geneLocationsIndexes)
        % The module is empty, so do not find its average
        continue
    end
    nonEmptyModules = vertcat(nonEmptyModules, k);
    for i = 1:nCellTypes
        % Iterate through cell types and find the average regulator
        % relevance for each
        DNaseRelAverage(:,k,i) = mean(DNaseRel(:,geneLocationsIndexes,i), 2)'; % Using 1 l1 penalty that is the average of all of the l1 penalties for the genes in the module
    end
end

weightTimesDNaseRelAverage = zeros(length(nonEmptyModules), size(DNaseRelAverage,1), size(DNaseRelAverage, 3));
for i = 1:nCellTypes
    % Iterate through cell types and multiply the absolute value of the 
    % weight by the meta-prior for each cell type
    weightTimesDNaseRelAverage(:,:,i) = abs(W(nonEmptyModules,:)) .* (DNaseRelAverage(:,nonEmptyModules, i)');
end
objectiveWithMetapriorVarsCellType = zeros(length(nonEmptyModules), size(DNaseRelAverage,1));
for ct = 1:length(cellTypeMap)
    % Iterate through all cell types and add the information for the prior
    % for the nearest cell type with DNase to the inner part of the
    % objective
    objectiveWithMetapriorVarsCellType = objectiveWithMetapriorVarsCellType + weightTimesDNaseRelAverage(:,:,cellTypeMap(ct)) - permute(log(DNaseRelAverage(:,nonEmptyModules,cellTypeMap(ct))), [2 1 3]);
end
objectiveWithMetapriorVars = sum(sum(objectiveWithMetapriorVarsCellType, 2)) + (e*sum(DNasePairWeights.^2));

geneGrad = {};
for k = 1:length(DNasePairFeats)
    % Initialize gene gradient for each feature
    geneGrad{k} = zeros(size(DNasePairFeats{1}));
end
for m = 1 : ngenes;
    % Compute the gradient with respect to the weights for each gene
    rel = zeros( size(DNasePairFeats{1},2), nCellTypes);
    for k = 1 : length(DNasePairFeats)
        % Iterate through pairwise features and add them appropriately
        % Assumes that pairwise features are 0 for DNases that are not
        % relevant to the gene
        for i = 1:nCellTypes
            % Iterate through cell types and find the relevance for each cell
            % type
            rel(:,i) = rel(:,i) + DNasePairFeats{k}(m,:,i)' * DNasePairWeights(k);
        end
    end
    z = computeSigmoidMat( rel );
    for k = 1:length(DNasePairFeats)
        % Iterate through features and compute the gradient for each
        geneGrad{k}(m,:,:) = (l1Lower-l1Upper)*(z .* (1-z).*permute(DNasePairFeats{k}(m,:,:), [2 3 1]));
    end
end

pgrad = 2 * e * DNasePairWeights;

for j = 1:size(u,2)
    % Iterate through modules and find the gradient for each module
    if isempty(find(ismember(nonEmptyModules, j)))
        % There are no genes in the current module, so continue
        continue
    end
    geneLocationIndexes = find(ismember(assign, j));
    if isempty(geneLocationIndexes)
        % The module is empty, so do not find the average of the gradient
        % corresponding to it
        continue
    end
    for k = 1 : length(DNasePairFeats)
        % Compute the gradient with respect to the weight for each feature
        % Compute the derivative of the prior, accounting for the Laplacian
        % normalization constant
        uTimesGradAverage = zeros(size(u,1), nCellTypes);
        gradAverage = zeros(size(u,1), nCellTypes);
        for i = 1:nCellTypes
            % Iterate through the cell types and find the term in the
            % gradient for each
            gradAverage(:,i) = mean(geneGrad{k}(geneLocationIndexes,:,i), 1);
            uTimesGradAverage(:,i) = sum(u(:,j).*gradAverage(:,i));
        end
        pgradkCellType = zeros(size(u,1), 1);
        for ct = 1:length(cellTypeMap)
            % Iterate through all cell types and add the gradient
            % information for the prior for the nearest cell type with
            % DNase to the inner part of the gradient
            pgradkCellType = pgradkCellType + uTimesGradAverage(:,cellTypeMap(ct)) - (gradAverage(:,cellTypeMap(ct))./permute(DNaseRelAverage(:,j,cellTypeMap(ct)), [1 3 2]));
        end
        pgrad(k) = pgrad(k) + sum(pgradkCellType);
    end
end
end


function [projectionOntoNegative, distNorm] = projectOntoNegative(featWeights)
% Computes the projection of feature weights onto the set of POSITIVE
% numbers
distSquared = 0;
for i = 1:length(featWeights)
    % Iterate through feature weights and set all negative feature weights
    % to 0
    if featWeights(i) < 0
        % The weight is negative, so set it to 0
        distSquared = distSquared + (abs(featWeights(i))^2);
        featWeights(i) = 0;
    end
end
projectionOntoNegative = featWeights;
distNorm = sqrt(distSquared);
end