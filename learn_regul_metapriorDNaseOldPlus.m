function DNasePairWeights = learn_regul_metapriorDNase (DNasePairFeats, cellTypeMap, DNasePairWeights, W, assign, nonEmptyModules, l1Upper, l1Lower, e)

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
%       nonEmptyModules: list of modules that contain at least 1 gene
%       l1Upper: the upper bound on the l1 penalty (Will be multiplied by 
%           (1-sigmoid))
%       l1Lower: the lower bound on the l1 penalty (Will be added to the L1 
%           penalty)
%       num_mod: is the number of modules
%       e : penalty for l2 regularization for DNase feature weights
% output:
%       DNasePairWeights: new values of meta-feature weights

fprintf('BOF::learn_regul_metaprior\n');

Wsum = sum(abs(W), 2);
u = zeros(size(W'));
for i = 1:size(u, 2)
    % Iterate through modules and normalize the weights in each module so
    % that the size of the weights overall does not affect the meta-prior
    % weights
    if Wsum(i) > 0
        % There is at least 1 TF with a non-zero weight
        u(:,i) = abs(W(i,:)')/Wsum(i);
    end
end


% f = @(DNasePairWeights)computeObjectiveWithMetapriorVars(DNasePairWeights, DNasePairFeats, cellTypeMap, u, assign, nonEmptyModules, l1Upper, l1Lower, e);
% options = optimset('Algorithm', 'interior-point', 'Hessian', 'lbfgs', 'SubproblemAlgorithm', 'cg', 'GradObj','on', 'TolFun', 0, 'UseParallel', 'always');
% A = -1 * eye(length(DNasePairFeats));
% b = zeros(length(DNasePairFeats), 1);
% lb = zeros(length(DNasePairFeats), 1) + .000001;
% DNasePairWeights = fmincon(f, DNasePairWeights, A, b, [], [], lb, Inf, [], options);
[DNasePairWeights,~,~] = minConF_PQNDNase(@computeObjectiveWithMetapriorVars, DNasePairWeights, @projectOntoPositive, DNasePairFeats, cellTypeMap, u, assign, nonEmptyModules, l1Upper, l1Lower, e);

fprintf('EOF::learn_marker_metaprior\n\n');
end


function [objectiveWithMetapriorVars, pgrad] = computeObjectiveWithMetapriorVars(DNasePairWeights, DNasePairFeats, cellTypeMap, u, assign, nonEmptyModules, l1Upper, l1Lower, e)
% Computes the objective given a set of feature weights
nCellTypes = size(DNasePairFeats{1}, 3);

DNaseRel = calculate_regul_relevanceDNase (l1Upper, l1Lower, DNasePairFeats, DNasePairWeights);
DNaseRelAverage = zeros(size(u,1), size(u,2), nCellTypes);
for k = 1:length(nonEmptyModules)
    % Compute the average relevance for each regulator across genes in the
    % module
    geneLocationsIndexes = find(ismember(assign, nonEmptyModules(k)));
    DNaseRelAverage(:,nonEmptyModules(k),:) = permute(mean(DNaseRel(:,geneLocationsIndexes,:), 2), [2 1 3]);
end

weightTimesDNaseRelAverage = zeros(length(nonEmptyModules), size(DNaseRelAverage,1), size(DNaseRelAverage, 3));
for i = 1:nCellTypes
    % Iterate through cell types and multiply the absolute value of the 
    % weight by the meta-prior for each cell type
    weightTimesDNaseRelAverage(:,:,i) = u(:,nonEmptyModules)' .* (DNaseRelAverage(:,nonEmptyModules, i)');
end

objectiveWithMetapriorVarsCellType = sum(weightTimesDNaseRelAverage(:,:,cellTypeMap) - permute(log(DNaseRelAverage(:,nonEmptyModules,cellTypeMap)), [2 1 3]),3);
objectiveWithMetapriorVars = sum(sum(objectiveWithMetapriorVarsCellType, 2)) + (e*sum(DNasePairWeights.^2));

geneGrad = {};
for k = 1:length(DNasePairFeats)
    % Initialize gene gradient for each feature
    geneGrad{k} = zeros(size(DNasePairFeats{1}));
end

rel = zeros(size(DNasePairFeats{1}));
for k = 1:length(DNasePairFeats)
      % Iterate through pairwise features and add them appropriately
      % Assumes that pairwise features are 0 for DNases that are not
      % relevant to a gene
      rel = rel + (DNasePairFeats{k} * DNasePairWeights(k));
end

z = computeSigmoidMat(rel);
for k = 1:length(DNasePairFeats)
    % Iterate through features and compute the gradient for each
    geneGrad{k} = (-l1Upper)*(z .* (1-z).*DNasePairFeats{k});
end

pgrad = 2 * e * DNasePairWeights;

for j = 1:size(u,2)
    % Iterate through modules and find the gradient for each module
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
        gradAverage = permute(mean(geneGrad{k}(geneLocationIndexes,:,:), 1), [2 3 1]);
        for i = 1:nCellTypes
            % Iterate through the cell types and find the term in the
            % gradient for each
            uTimesGradAverage(:,i) = u(:,j).*gradAverage(:,i);
        end
        pgradkCellType = sum(uTimesGradAverage(:,cellTypeMap) - (gradAverage(:,cellTypeMap)./permute(DNaseRelAverage(:,j,cellTypeMap), [1 3 2])), 2);
        pgrad(k) = pgrad(k) + sum(pgradkCellType);
    end
end
end


function [projectionOntoPositive, distNorm] = projectOntoPositive(featWeights)
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
projectionOntoPositive = featWeights;
distNorm = sqrt(distSquared);
end