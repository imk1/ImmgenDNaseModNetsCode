function W = learn_weightsDNase (tf_expression, modmean, assign, DNaseRel, cellTypeMap, d)

% Code written by Su-In Lee
% Modified by Irene Kaplow

% This function learns the weigths for the regulators

% input:
%       tf_expression:  tf expression, size rxn where r is the number of
%           regulators, each TF's expression has been normalized (if necessary)
%       modmean: the average expression of the genes in each module for
%           each cell type -- matrix is num_mods x n, where num_mods is the
%           number of modules and n is the number of cell types
%       assign: module assignments
%       DNaseRel: the relevance of each regulator for each gene in each
%           cell type -- matrix that is r x g x i, where r is the number of
%           regulators, g is the number of genes, and i is the number of cell
%           types with DNase
%       cellTypeMap: maps each cell type to the closest cell type with
%           near-by DNase -- vector of length n, where n is the number of 
%           cell types and where each entry has the index of the closest
%           cell type with DNase
%       d: L2 parameter for feature weights
% output:
%       W: a matrix of size num_mods x r, that contains the 
%           coefficients for each module, regulator pair

fprintf('BOF::learn_weights\n');

nfeats = size(tf_expression,1);
[ncls, ~] = size(modmean);

feat = tf_expression;
featT = feat';

W = zeros(ncls, nfeats);

for k = 1 : ncls
    % Iterate through modules

   if( mod(k,40)==0 )
      % Print iteration information
      fprintf('k = %g / ncls = %g\n',k,ncls);
   end
   
   geneLocationsIndexes = find(ismember(assign, k));
   if isempty(geneLocationsIndexes)
       % The module is empty, so do not find its weights
        continue
   end
   l1 = mean(DNaseRel(:,geneLocationsIndexes,cellTypeMap), 2);
   currentModMean = modmean(k, :);
   W(k,:) = elasticNetADMM(featT, currentModMean', permute(l1, [1 3 2]), 1, 1, d);

end

fprintf('EOF::learn_weights\n\n');

