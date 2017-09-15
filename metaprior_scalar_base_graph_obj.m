function [o,g] = metaprior_scalar_base_graph_obj(A,F,W, L,c1,c2);
% function [o,grad] =metaprior_scalar_obj(A,F,W,c1);
U = abs(W);
[r,m] = size(U);
[r,p] = size(F);
sc = (A(1));
bs = A(2);
A = A(3:end);
A = reshape(A,p,m);
rel = sigmoid(F*A);

o = sum(sum(sc.*(bs+rel).*U)) - sum(sum(log(sc.*(bs+rel)))) + ...
       c1*trace(A*L*A') +  sum(c2.*sum(A.^2,2)) ;
g = F'*(sc.*U.*rel.*(1-rel)-(1./(bs+rel)).*rel.*(1-rel)) + 2*diag(c2)*A + 2*c1*A*L;
g2 = sum(sum((bs+rel).*U)) - (1/sc)*(r*m);
g3 = sc.*sum(sum(U)) - sum(sum(1./(bs+rel)));
g = [g2;g3;g(:)];



