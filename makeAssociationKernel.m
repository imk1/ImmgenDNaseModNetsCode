function G = makeAssociationKernel(X, K)
%function G =  makeAssociationKernel(X, K)
%
%  input: 
%   X is a MxN matrix of data
%   Columns of X are the datapoints, i.e. there are N genes in X and each gene has
%   measurements for M features.
%   K is the number of neighbours that are kept for each gene (K should be
%   betwee 1 and N - 1);
%  
%  Returns a NxN sparse matrix G
%  G(I,J) = 0 if J is not within the K closest neighbours
%  otherwise, G(I,J) = the Pearson correlation of X(:,I) and X(:,J)

[M,N] = size(X);
G = sparse([],[],[],N,N,N*(K+1));
XX = X - repmat(mean(X, 1), M, 1);
XX = XX ./ repmat(sqrt(M)*std(XX, 1, 1), M, 1);
XXt = XX';
fprintf('makeKernel called with %d datapoints.  One dot / 100 datapoints\n', N);
for ii = 1:N
    value = XXt(ii,:) * XX;
    [sortedValues, indSort] = sort(-1*value);
    G(indSort(1:K+1), ii) = -1*sortedValues(1:K+1)';
%    G(ii, indSort(1:K+1)) = sortedValues(1:K+1);
    if mod(ii, 100) == 0
        fprintf('.');
        if mod(ii, 8000) == 0
            fprintf('\n');
        end
    end
end

G = G - diag(diag(G));
G = max(G, G');
