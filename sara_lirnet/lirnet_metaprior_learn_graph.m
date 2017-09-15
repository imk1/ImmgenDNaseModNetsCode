function [W,A,obj] = lirnet_metaprior_learn_graph(modmeans,genotype,features,Ln,Winit,lambda1,lambda2,lambda3,iter);
% [W,A,obj,scalar,base] =
% lirnet_metaprior_learn_graph(modmeans,genotype,features,Ln,Winit,lambda1,lambda2,lambda3);

[nsamples,nregs] = size(genotype);
num_modules = size(modmeans,2);
nfeatures = size(features,2);
rho = 1;
alpha = 1;
verbose = 0;
W = Winit;
A = rand(nfeatures,num_modules);
mybase = 0.001;


for ii = 1:iter
    o = zeros(num_modules,1); 

    [AA,A] = solve_metaprior_scalar_graph(features,W,Ln,lambda2,lambda3);
    rel = sigmoid(features*A);
    
    
    for jj = 1:num_modules
        lambdas = mybase+lambda1*(rel(:,jj));
        z = lassoADMMSolver(genotype, modmeans(:,jj), lambdas,0,rho,alpha,verbose);
        Wnew(:,jj) = z;
        o(jj) = sum((genotype*z-modmeans(:,jj)).^2) + sum(sum(abs(z).*lambdas));
    end
    W = Wnew;

    obj(ii) = sum(o);
end


        
        
        
    



