function zList = runElasticNetADMM(AFileName, bFileName, rho, alpha, a, zListOutputFileName, historyListOutputFileName, residualListOutputFileName, predictionListOutputFileName)
% Runs the ADMM version of the elastic net with multiple lambda values for
% multiple developmental stages

lambdaList = [0; .0000001; .000001; .00001; .0001; .001];
A = importdata(AFileName);
bMat = importdata(bFileName);

% Order of lists: first lambda for first stage, second lambda for first
% stage, ..., last lambda for first stage, first lambda for second stage,
% ... last lambda for last stage
zList = zeros(size(A,2), size(bMat, 2)*length(lambdaList));
historyList = {};
residualList = zeros(size(A,1), size(bMat, 2)*length(lambdaList));
predictionList = zeros(size(A,1), size(bMat, 2)*length(lambdaList));

% Make sure number input is doubles and not strings
rho = str2num(rho);
alpha = str2num(alpha);
a = str2num(a);

entryCount = 0;
for i = 1:size(bMat, 2)
    % Iterate through developmental stages
    b = bMat(:, i);
    for j = 1:length(lambdaList)
        % Iterate through lambdas
        entryCount = entryCount + 1;
        lambda = lambdaList(j);
        [z, history] = elasticNetADMM(A, b, lambda, rho, alpha, a);
        zList(:,entryCount) = z;
        historyList{entryCount} = history;
        clear history
        predictions = A*z;
        clear z
        residuals = predictions - b;
        residualList(:,entryCount) = residuals;
        clear residuals
        predictionList(:, entryCount) = predictions;
        clear predictions
    end
    clear b
end

save (zListOutputFileName, 'zList');
save (historyListOutputFileName, 'historyList');
clear historyList
save (residualListOutputFileName, 'residualList');
save (predictionListOutputFileName, 'predictionList');