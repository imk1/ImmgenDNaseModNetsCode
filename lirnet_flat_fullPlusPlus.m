function params = lirnet_flat_fullPlusPlus(expression,tf_expression,l1,l2)

% Code written by Sara Mostafavi
% Modified by Irene Kaplow

% input:
%       expression: the expression file, size gxn where g is the number of genes
%               and n is the number of patients
%       tf_expression:  tf expression, size rxn wherer is the number of regulators
%       l1 : penalty factor for l1 regularization
%       l2 : penalty factor for l2 regularization
%       num_mod: is the number of modules
% output:
%       params is a structure with several relevant fields:
%           errors: is the sum squared error
%           weights: is a matrix of size r x num_mods , coefficients for each module
%           assign: has d columns, each column is the assigment of genes into cluster
%               for the given iteration (d iteration in total) (note: do not use more
%               than 2 iterations if the penalty for l1 regularization is stringent).

genotype = tf_expression;
[nreg,narrays] = size(genotype);
[ngenes,nsamples] = size(expression);


[num_genes,num_arrays] = size(expression);
pp = randperm(narrays);
reclust =3;

G = makeAssociationKernel(expression', 100);
[idx,netsim,dpsim,expref]=apcluster(G,0);
indexes = unique(idx);
assign = zeros(length(idx), 1);
for i = 1:length(idx)
    index = find(indexes == idx(i));
    assign(i) = index;
end
%assign = kmeans_clustering(expression,num_mod,10);

W = zeros(nreg,max(assign));

for rr = 1:reclust;
    % Learn the regression for each module for each re-assignment of genes
    % to clusters
    
    dist = Inf(num_genes,max(assign));
    modmean = Inf(max(assign),nsamples);
    
    for ii = 1:max(assign)
        % Iterate through modules and learn the regression for each module
        
        g = find(assign==ii);
        if length(g)==0
            % Do not do anything for an empyt module
            continue
        elseif length(g)>1
            % Find the mean expressions across the genes in the module
            modmean(ii,:) = mean(expression(g,:));
        else
            modmean(ii,:) = expression(g,:);
        end
        
        y = modmean(ii,:)';
        
        [W(:,ii), history] = elasticNetADMM(genotype', y, l1, 1, 1, l2);
        pred = genotype'*W(:,ii); % predictions
        params.errors(rr,ii) = sum ( sum ( (expression(g,:) -  repmat(pred',length(g),1)).^2)); %sse
        m = mean(expression(g,:),2);
        
        dist(:,ii) = sum ( (expression - repmat(pred',num_genes,1)).^2,2);
    
    end
    
    [junk,new_assign] = min(dist,[],2);
    nmoved = sum(new_assign~=assign);
    fprintf('genes moved %d\n',nmoved);
    
    assign = new_assign;
    params.assign(:,rr) = assign;
    params.weights{rr} = W;
    if nmoved < 0.05*(num_genes);
        break
    end
    
end












