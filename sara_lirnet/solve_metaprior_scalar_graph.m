function [Anew,a,z] = solve_metaprior_scalar_graph(F,W,L,c1,c2);

U = abs(W);
[r,m] = size(U);
[r,p] = size(F);
A = zeros(p,m);
A = [10;0.1;A(:)];
options.verbose = 0;
options.display = 0;
options.maxIter = 100;
options.method = 'lbfgs';
funObj = @(w)metaprior_scalar_base_graph_obj(w,F,W,L,c1,c2);
funProj = @(w)project_positive(w);
[Anew ,z]= minConF_PQN(funObj,A,funProj,options);
a = reshape(Anew(3:end),p,m);




end



function p = project_positive(r);

p = r;
p = p(:);
if(p(1))<0
    p(1) = 0;
 end
 p3 = p(3:end);
 ix = find(p3>0);
 p3(ix) = 0;
 p = [p(1:2);p3(:)];
 end
