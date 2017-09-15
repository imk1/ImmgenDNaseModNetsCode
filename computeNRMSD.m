function NRMSDs = computeNRMSD(data, featureData, bMat)
% Computes the NRMSD for a multiple regressions, where each row in data had
% its own regression

preds = zeros(size(data));
for i = 1:size(data, 1)
    % Iterate through genes and find the predicted value for each gene
    weights = bMat(i, :);
    for j = 1:size(data, 2)
        % Iterate through cell types and find the predicted value for each
        % cell type for the current gene
        %preds(i, j) = weights(1) + sum(weights(2:length(weights)) .* featureData(:, j)');
        preds(i, j) = sum(weights .* featureData(:, j)');
    end
end

resids = zeros(size(data, 1));
for i = 1:size(data, 1)
    % Find the residual for the regression for each gene
    correctData = data(i,:);
    resids(i) = sqrt(sum(((correctData - preds(i, :)).^2)));
end

NRMSDs = zeros(size(data, 1) , 1);
sqrtn = sqrt(size(data, 2));
for i = 1:size(data, 1)
    % Find the NRMSD for the regression for each gene
    NRMSDs(i) = (resids(i) / sqrtn) / (max(data(i,:)) - min(data(i,:)));
end