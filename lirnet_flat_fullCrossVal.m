function params = lirnet_flat_fullCrossVal(expressionFull,tf_expressionFull,assign, num_mod,l1,l2, validationIndexes, reclust)

% Code written by Sara Mostafavi
% Modified by Irene Kaplow

% input:
%       expressionFull: the expression matrix, size gxn where g is the 
%           number of genes and n is the number of cell types
%       tf_expressionFull:  tf expression, size rxn wherer is the number of 
%           regulators
%       assign: initial assignments of genes to modules
%       num_mod: is the number of modules
%       l1 : penalty factor for l1 regularization
%       l2 : penalty factor for l2 regularization
%       validationIndexes: the indexes of the cell types that will be used
%           for validation
%       reclust: the number of times that re-clustering will occur
% output:
%       params is a structure with several relevant fields:
%           errors: is the sum squared error for the training set
%           weights: is a matrix of size r x num_mods , coefficients for 
%               each module
%           assign: has d columns, each column is the assigment of genes 
%               into cluster for the given iteration (d iteration in total) 
%               (note: do not use more than 2 iterations if the penalty for 
%               l1 regularization is stringent).
%           errorsVal: is the sum squared error for the validation set

trainingIndexes = setdiff(1:size(expressionFull, 2), validationIndexes);
expression = expressionFull(:, trainingIndexes);
tf_expression = tf_expressionFull(:, trainingIndexes);
expressionVal = expressionFull(:, validationIndexes);
tf_expressionVal = tf_expressionFull(:, validationIndexes);

[nreg,~] = size(tf_expression);
[num_genes,nsamples] = size(expression);

W = zeros(nreg,num_mod);

for rr = 1:reclust;
    % Learn the regression for each module for each re-assignment of genes
    % to clusters
    
    dist = Inf(num_genes,num_mod);
    modmean = Inf(num_mod,nsamples);
    
    for ii = 1:num_mod
        % Iterate through modules and learn the regression for each module
        
        g = find(assign==ii);
        if length(g)==0
            % Do not do anything for an empty module
            continue
        elseif length(g)>1
            % Find the mean expressions across the genes in the module
            modmean(ii,:) = mean(expression(g,:));
        else
            modmean(ii,:) = expression(g,:);
        end
        
        y = modmean(ii,:)';
        
        [W(:,ii), ~] = elasticNetADMM(tf_expression', y, l1, 1, 1, l2);
        pred = tf_expression'*W(:,ii); % predictions
        params.errors(rr,ii) = sum ( sum ( (expression(g,:) -  repmat(pred',length(g),1)).^2)); %sse
        
        predVal = tf_expressionVal'*W(:,ii);
        params.errorsVal(rr,ii) = sum ( sum ( (expressionVal(g,:) -  repmat(predVal',length(g),1)).^2)); %sse
        
        dist(:,ii) = sum ( (expression - repmat(pred',num_genes,1)).^2,2);
    
    end
    
    [~,new_assign] = min(dist,[],2);
    nmoved = sum(new_assign~=assign);
    fprintf('genes moved %d\n',nmoved);
    
    assign = new_assign;
    params.assign(:,rr) = assign;
    params.weights{rr} = W;
    if nmoved < 0.05*(num_genes);
        break
    end
    
end












