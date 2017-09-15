function X = standardize(X)

% Code written by Sara Mostafavi
% Modified by Irene Kaplow

% Normalizes each row of X

[n p] = size(X);

X = X - nanmean(X,2)*ones(1,p);

X(isnan(X)) = 0;

X = X./sqrt(sum(X.^2,2)*ones(1,p));