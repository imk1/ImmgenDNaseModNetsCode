function sigVal = computeSigmoidMat(z)
% This function computes the value of the sigmoid of all of the numbers in 
% an n x m x p matrix, where n, m, and p are the dimensions of z.
%
% Input:
%   z: The n x m x p matrix of values whose 
%   sigmoids need to be found
%
% Output:
%   sigVal: The value of the sigmoid

expNegz = exp(-z);
sigVal = 1 ./ (1 + expNegz);