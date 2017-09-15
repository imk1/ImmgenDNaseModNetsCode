function DNasePairWeights = learn_regul_metapriorDNase (DNasePairFeats, DNasePairWeights, W, assign, l1Upper, l1Lower, e, PWMScores)

% Code written by Su-In Lee
% Modified by Irene Kaplow

% Calculates the meta-priors for the DNase features

fprintf('BOF::learn_regul_metaprior\n');

opts = optimset('Algorithm', 'interior-point');
problem = struct();
problem.objective = @computeObjectiveWithMetapriorVars;
problem.x0 = DNasePairWeights;
problem.Aineq = 1;
problem.Bineq = zeros(length(DNasePairWeights), 1);
problem.solver = 'fmincon';
problem.options = opts;
DNasePairWeights = fmincon(problem); % Minimize the objective subject to each weight being <= 0

fprintf('EOF::learn_marker_metaprior\n\n');
end

function objectiveWithMetapriorVars = computeObjectiveWithMetapriorVars(featWeights)
% Computes the objective given a set of feature weights
DNaseRel = calculate_regul_relevanceDNase (l1Upper, l1Lower, DNasePairFeats, featWeights, PWMScores);
DNaseRelAverage = zeros(size(u));
for k = 1:size(u, 2)
    % Compute the average relevance for each regulator across modules
    geneLocationsIndexes = find(ismember(assign, k));
    DNaseRelAverage(:,k) = mean(DNaseRel(:,geneLocationsIndexes), 2)'; % Using 1 l1 penalty that is the average of all of the l1 penalties for the genes in the module
end
objectiveWithMetapriorVars = sum(sum( abs(W) .* DNaseRelAverage' - log(DNaseRelAverage'), 2)) + (e*sum(featWeights.^2));
end