function idx = hierarchicalClusterDivisive(dataMat, clustStopSize, iterStop)
% Takes a matrix that is examples by features and does hierarchical
% divisive clustering, stopping when clusters are less than or equal to a
% specified size; the return value is the list of cluster numbers for each
% example

idx = zeros(size(dataMat, 1), 1);
dataQueue = {};
dataQueue{1} = 1:1:size(dataMat, 1);
clustCount = 0;
iterCount = 0;

while ~isempty(dataQueue)
    % Keep clustering until all examples have been assigned a final cluster
    iterCount = iterCount + 1
    dataCurrentLocations = dataQueue{1};
    dataCurrent = dataMat(dataCurrentLocations, :);
    [~,C] = kmeansPlusPlus(dataCurrent',2);
    idxTemp = kmeans(dataCurrent, 2, 'start', C');
    %idxTemp = kmeans(dataCurrent, 2, 'replicates', 20);
    dataQueue = dataQueue(2:length(dataQueue));
    for i = 1:2
        % Iterate through current clusters to see if any of them are final
        % clusters
        if length(find(idxTemp == i)) <= clustStopSize
            % Assign each example in the first cluster to its final cluster
            % number
            clustCount = clustCount + 1;
            for j = 1:length(idxTemp)
                % Iterate through examples that are assigned to their final
                % cluster
                if idxTemp(j) == i
                    % Assign the example to its final cluster number
                    trueLocation = dataCurrentLocations(j);
                    idx(trueLocation) = clustCount;
                end
            end
        else
            dataQueue{length(dataQueue) + 1} = dataCurrentLocations(idxTemp == i);
        end
    end
    if iterCount == iterStop
        % Maximum number of iterations has happened
        break
    end
end

if iterCount == iterStop
    % Stopped because the maximum number of iterations had occurred
    for i = 1:length(dataQueue)
        % Iterate through clusters and make cluster number assignments
        dataCurrentLocations = dataQueue{i};
        clustCount = clustCount + 1;
        for j = 1:length(dataCurrentLocations)
            % Assign each example to its final cluster number
            trueLocation = dataCurrentLocations(j);
            idx(trueLocation) = clustCount;
        end
    end
end
        