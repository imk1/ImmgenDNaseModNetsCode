function centroidMat = getModuleCentroids(exprData, assign)
% Gets the centroids of modules

centroidMat = zeros(max(assign), size(exprData, 2));

for i = 1:max(assign)
    % Iterate through the modules and compute the centroid for each
    genesInModule = find(assign == i);
    if length(genesInModule) > 0
        % The module is not empty, so compute the centroid
        exprModule = exprData(genesInModule,:);
        centroidMat(i,:) = mean(exprModule, 1);
    end
end