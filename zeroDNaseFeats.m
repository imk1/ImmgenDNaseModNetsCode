function DNasePairFeats = zeroDNaseFeats(DNaseFeatFileNames, binaryFeatIndexes)
% This function loads DNase feature matrices into Matlab, sets all of the
% -1's to 0's, sets all of the DNase, gene feature pairs to 0 if any of
% them are -1, and puts all of the matrices into a cell array

for i = 1:length(DNaseFeatFileNames)
    % Iterate through feature matrices and import all of them into Matlab
    % Iterate through feature matrices and modify -1's and 0's
    % appropriately
    i
    DNasePairFeats{i} = importdata(DNaseFeatFileNames{i});
    for j = 1:size(DNasePairFeats{i}, 1)
        % Iterate through DNase sites
        for k = 1:size(DNasePairFeats{i}, 2)
            % Iterate through genes
            if DNasePairFeats{i}(j,k) == -1
                % DNase feature is not available, either because the DNase
                % is not sufficiently close to the gene's TSS or because
                % the signal was not significant, so change the -1 to a 0
                % In addition, change (j,k)th entry in the earlier feature
                % matrices to a 0
                % Also, change the (j,k)th entry in the later feature
                % matrices that are not binary features to -1
                DNasePairFeats{i}(j,k) = 0;
                for l = 1:i-1
                    % Iterate through earlier feature matrices and change
                    % the (j,k)th entries to 0
                    DNasePairFeats{l}(j,k) = 0;
                end
                for l = i+1:length(DNasePairFeats)
                    % Iterate through later feature matrices and change the
                    % (j,k)th entries to -1 in matrices for non-binary
                    % features
                    if isempty(find(binaryFeatIndexes == l))
                        % The feature matrix is not for a binary feature,
                        % so change the (j,k)th entry to -1
                        DNasePairFeats{l}(j,k) = -1;
                    end
                end
            end
        end
    end
    DNasePairFeats{i} = sparse(DNasePairFeats{i});
end

save ../MouseDNase/DNasePairFeatsSparse DNasePairFeats