function [X] = whiten(X)

% Code written by Patrick J. Mineault and found on 
% http://xcorr.net/2011/05/27/whiten-a-matrix-matlab-code/
% Modified by Irene Kaplow

% This function transforms the matrix X into a matrix that has identity
% covariance matrix

X = bsxfun(@minus, X, mean(X));
A = X'*X;
[V,D] = eig(A);
X = X*V*diag(1./(diag(D)).^(1/2))*V';