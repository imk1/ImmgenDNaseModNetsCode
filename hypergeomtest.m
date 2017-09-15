
function pp = hypergeomtest(kk, nn, CC, GG)
%HYPERGEOM - hypergeometric p-value
%
%     P = HYPERGEOM(K, N, C, G)
%
%  Probability of seeing at least K members of a set of size N in a
%  cluster of size C of G total entities.
%


%  C(N,K) = LOG N choose K
C = inline('gammaln(N+1) - gammaln(K+1) - gammaln(N-K+1)', ...
          'N', 'K');

tmp = [];

upper = min(nn, CC);
for ii = kk:upper
 tmp(ii) = C(CC,ii) + C(GG-CC, nn-ii) - C(GG, nn);
end

pp = sum(exp(tmp(upper:-1:kk)));