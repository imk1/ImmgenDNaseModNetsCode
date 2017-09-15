function DNaseRel = calculate_regul_relevanceDNase (l1Upper, l1Lower, DNasePairFeats, DNasePairWeights)

% Code written by Su-In Lee
% Modified by Irene Kaplow

% This function computes the l1 penalty based on the metaprior weights for
% the DNase features

% input:
%       l1Upper: the upper bound on the l1 penalty (Will be multiplied by 
%           (1-sigmoid))
%       l1Lower: the lower bound on the l1 penalty (Will be added to the L1 
%           penalty)
%       DNasePairFeats: pairwise metafeatures for DNase sites -- cell array
%           where each cell has a different pairwise feature for a gene and 
%           a TF, each cell's entry is g x r x i, where g is the number of
%           genes, r is the number of regulators (assumes that all 
%           features relating to a DNase site that are not near a gene have 
%           value 0 and that PWM scores have already been incorporated into
%           features), and i is the number of cell types with DNase
%       DNasePairWeights: initial values of meta-feature weights
% output:
%       DNaseRel: the relevance of each regulator for each gene in each
%           cell type -- matrix that is r x g x i, where r is the number of
%           regulators, g is the number of genes, and i is the number of cell
%           types with DNase

ngenes = size( DNasePairFeats{1},1 );
nreguls = size(DNasePairFeats{1}, 2);
nCellTypes = size(DNasePairFeats{1}, 3);
nFeatTypes = length(DNasePairFeats);

rel = zeros(nreguls, ngenes, nCellTypes);

for k = 1:nFeatTypes
      % Iterate through pairwise features and add them appropriately
      % Assumes that pairwise features are 0 for DNases that are not
      % relevant to a gene
      for j = 1:nreguls
          % Iterate through TFs and add the relevance for each
          rel(j,:,:) = rel(j,:,:) + (permute(DNasePairFeats{k}(:,j,:), [2 1 3]) * DNasePairWeights(k,j));
      end
end
DNaseRel = l1Upper * (ones(size(rel)) - computeSigmoidMat( rel )) + l1Lower; % Number of regulators by number of genes by number of cell types with DNase