function W = learn_weightsDNase (tf_expression, modmean, assign, DNaseRel)

% Code written by Su-In Lee
% Modified by Irene Kaplow

% This function learns the weigths for the regulators

fprintf('BOF::learn_weights\n');

nfeats = size(tf_expression,1);
[ncls, nsamples] = size(modmean);

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
   l1 = mean(DNaseRel(:,geneLocationsIndexes), 2); % Using 1 l1 penalty that is the average of all of the l1 penalties for the genes in the module
   currentModMean = modmean(k, :);
   W(k,:) = lassoADMM(featT, currentModMean', l1, 1, 1);

end

fprintf('EOF::learn_weights\n\n');

