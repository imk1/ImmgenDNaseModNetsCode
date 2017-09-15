function params = lirnet_metaprior_full_graph_wboot(expression,genotype,features,assign,l1param,l2graph,l2met,sgraph,gamma);
% params =
% lirnet_metaprior_full_graph_wboot(expression,genotype,features,assign,l1param,l2graph,l2met,sgraph,gamma);
% input:
%       expression: the expression file, size gxn where g is the number of genes
%               and n is the number of patients
%       genotype:  the covariates, size rxn wherer is the number of regulators
%       features: is the features for regulators , a matrix of rxz where r is the
%       number of regulators and z is the number of features (for now I artificially
%                  makes this into 'module-specific' features
%       c0 and c1 : penalty factor for l1 regularization (c0 is the minimum and c1 is
%       the maximum)
%       ereg : penalty factor for l2 regularization for metaprior
%       dreg : penalty factor for l2 regularization for weights (W)
%       num_mod: is the number of modules
% output:
%       params is a structure with several relevant fields:
%           errors is the sum squared error on test (one for each cross-vlidation
%           iteration)
%           meanvar is the null deviation
%           errorT is the training error
%           weights is a matrix of size r x num_mods , coefficients for each module
%           assign has d columns, each column is the assigment of genes into cluster
%           for the given iteration (d iteration in total) (note: do not use more
%           than 2 iterations if the penalty for l1 regularization is stringent).
%           metaprior: is the coefficients for metafeatures

[nreg,narrays] = size(genotype);
[ngenes,nsamples] = size(expression);
reclust = 2;


num_mod = length(unique(assign));
params.errors = NaN(reclust,num_mod);

% get training indices
[num_genes,num_arrays] = size(expression);



for rr = 1:reclust;
    
    dist = Inf(num_genes,num_mod);
    modmean = Inf(num_mod,nsamples);
    W = zeros(nreg,num_mod);
    u = unique(assign);
    
    for ii = 1:num_mod
        
        g = find(assign==u(ii));
        if isempty(g)
            continue
        elseif length(g)>1
            modmean(ii,:) = mean(expression(g,:));
        else
            modmean(ii,:) = expression(g,:);
        end
        if rr==1
            Winit(:,ii) = lassoADMMSolver(genotype',modmean(ii,:)',gamma,0,1,1,0);
        end
       
    end
    CC = corr(modmean');
    CC(isnan(CC)) = 0;
    CC(CC<0) = 0;
    Cexp = CC.^sgraph;
    Ln = speye(num_mod)-normalizeKernel(Cexp);
   
    [W,MP] = lirnet_metaprior_learn_graph(modmean',genotype',features,Ln,Winit,l1param,l2graph,l2met,3);
    params.weights{rr} = W;
    params.metaprior{rr} = MP;
    
    for ii = 1:num_mod
        g = find(assign==ii);
        if isempty(g)
            continue
        end
        pred = genotype'*W(:,ii); 
        params.errors(rr,ii) = sum ( sum ( (expression(g,:) -  repmat(pred',length(g),1)).^2)); %sse
        dist(:,ii) = sum ( (expression - repmat(pred',num_genes,1)).^2,2);
    end
    
    [junk,new_assign] = min(dist,[],2);
    nmoved = sum(new_assign(:)~=assign(:));
    fprintf('genes moved %d\n',nmoved);
    
    assign = new_assign;
    params.assign(:,rr) = assign;
    params.weights{rr} = W;
    if nmoved < 0.05*(num_genes);
        break
    end
    
end




        
        
        
    



