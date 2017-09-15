function m = mcl(a, r, l, maxiters)
% function m = mcl(a, r, l, maxiters)
%
% run mcl algoirthm, on graph with weighted adj matrix a,
% multiplication exponent r, with 2 for the first l iters
% TAKEN FROM: http://www.cs.yale.edu/homes/spielman/462/2006/matlab/mcl.m,
% written by Dan Spielman

if (nargin < 2),
    r = 2;
end

if (nargin < 3),
    l = 0;
end

if (nargin < 4),
    maxiters = 50;
end

m = (a + diag(sum(a)));

i = 1;
while and(i < maxiters, sum(sum(abs(m - m^2))) > 10^(-6)),
    
  d = diag(1./sum(m));
  m = m * d;
  m = m ^ 2;
  
  if (i <= l)
    m = m .^ 2;
  else
    m = m .^ r;
  end
  
  i = i + 1;
  
end