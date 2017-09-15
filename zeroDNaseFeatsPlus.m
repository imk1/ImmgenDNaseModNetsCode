function DNasePairFeats = zeroDNaseFeatsPlus(DNaseFeatFileNames, binaryFeatIndexes)
% This function loads DNase feature matrices into Matlab, sets all of the 
% DNase, gene feature pairs to 0 if any of those that are for non-binary 
% features are 0, and puts all of the matrices into a cell array
% ASSUMES THAT, FOR EACH BINARY FEATURE, IF THE BINARY FEATURE IS NOT 
% PRESENT, SOME OTHER FEATURE WILL NOT BE PRESENT

for i = 1:length(DNaseFeatFileNames)
    % Iterate through feature matrices and import all of them into Matlab
    i
    DNasePairFeats{i} = importdata(DNaseFeatFileNames{i});
    DNasePairFeats{i} = sparse(DNasePairFeats{i});
end
for i = 1:length(DNasePairFeats)
    % Iterate through feature matrices and modify -1's and 0's
    % appropriately
    i
    if ~isempty(find(binaryFeatIndexes == i))
        % The feature matrix is for a binary feature, so do not use its
        % zero-values to determine when data is not present
        continue
    end
    for j = 1:size(DNasePairFeats{i}, 1)
        % Iterate through DNase sites
        for k = 1:size(DNasePairFeats{i}, 2)
            % Iterate through genes
            if DNasePairFeats{i}(j,k) == 0
                % DNase feature is not available, either because the DNase
                % is not sufficiently close to the gene's TSS or because
                % the signal was not significant, so change the -1 to a 0
                % In addition, change (j,k)th entry in the earlier feature
                % matrices to a 0
                % Also, change the (j,k)th entry in the later feature
                % matrices that are not binary features to -1
                for l = 1:i-1
                    % Iterate through earlier feature matrices and change
                    % the (j,k)th entries to 0
                    DNasePairFeats{l}(j,k) = 0;
                end
                for l = i+1:length(DNasePairFeats)
                    % Iterate through later feature matrices and change the 
                    % (j,k)th entries to 0
                    DNasePairFeats{l}(j,k) = 0;
                end
            end
        end
    end
    
end

save ../MouseDNase/DNasePairFeatsSparse DNasePairFeats