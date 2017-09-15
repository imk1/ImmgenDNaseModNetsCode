function params = lirnetDNase(expression, tf_expression, assign, DNasePairFeats, cellTypeMap, DNaseRel, num_mod, maxiters, reclust, l1Bounds, e)

% Based on code written by Sara Mostafavi and Su-In Lee
% Code written (mostly) by Irene Kaplow

% This code runs Lirnet with DNase-related priors that are gene-specific

% input:
%       expression: the expression file, size gxn where g is the number of genes
%               and n is the number of cell types, each gene's expression has
%               been normalized (if necessary)
%       tf_expression:  tf expression, size rxn where r is the number of
%           regulators, each TF's expression has been normalized (if necessary)
%       assign: initial module assignments, leave it empty for k-means
%           clustering
%       DNasePairFeats: pairwise metafeatures for DNase sites -- cell array
%           where each cell has a different pairwise feature for a gene and 
%           a TF, each cell's entry is g x r x i, where g is the number of
%           genes, r is the number of regulators (assumes that all 
%           features relating to a DNase site that are not near a gene have 
%           value 0 and that PWM scores have already been incorporated into
%           features), and i is the number of cell types with DNase
%       cellTypeMap: maps each cell type to the closest cell type with
%           near-by DNase -- vector of length n, where n is the number of 
%           cell types and where each entry has the index of the closest
%           cell type with DNase
%       DNaseRel: the initial values of relevance of each regulator for 
%           each gene in each cell type -- matrix that is r x g x n, where 
%           r is the number of regulators, g is the number of genes, and n 
%           is the number of cell types with DNase, leave empty for no
%           initialization
%       num_mod: is the number of modules
%       maxiters: the maximum number of iterations for the metaprior loop
%       reclust: number of times to move genes to better modules and repeat
%           the Lirnet process
%       l1Bounds: an array with 2 entries: the upper bound on the l1
%           penalty followed by the lower bound on the l1 penalty
%       e : penalty for l2 regularization for DNase feature weights
% output:
%       params is a structure with several relevant fields:
%           W0: matrix of size num_mods x r, where are is the number of
%               regulators, with the initial learned weights for each
%               regulator
%           errors: is the sum squared error
%           assign: has reclust columns, each column is the assigment of genes into cluster
%               for the given iteration (reclust iterations in total) (note: do not use more
%               than 2 iterations if the penalty for l1 regularization is stringent).
%           weights: cell array in which each entry is a matrix of size 
%               num_mods x r, coefficients for each module, where the number of
%               entires is equal to the number of re-clustering events
%           DNaseRel: the relevance of each regulator for each gene in each
%               cell type -- matrix that is r x g x n, where r is the number of
%               regulators, g is the number of genes, and n is the number of cell
%               types with DNase
%           DNasePairWeights: cell array in which each entry is a matrix of size hp, where hp is the number of
%               DNase pairwise features, that contains the learned
%               metapriors for the DNase pairwise features, and there is an
%               entry for each re-clustering iteration

% TF FEATURES NOT INCLUDED, ADD THEM LATER
% GENERAL L1 for TFS WITH NO PWM INFORMATION NOT INCLUDED, ADD IT LATER
[nregs,nsamples] = size(tf_expression);
num_genes = size(expression, 1);

if isempty(assign)
    % No initial modules have been given, so use k-means clustering to
    % assign initial modules
    if num_mod < size(expression, 1)
        % Learning module networks and not just a regression, so cluster the
        % genes
        [~, C] = kmeansPlusPlus(expression', num_mod);
        assign = kmeans(expression, num_mod, 'start', C');
    end
end

l1Upper = l1Bounds(1); % Will be multiplied by (1-sigmoid)
l1Lower = l1Bounds(2); % Will be added to L1 penalty

DNasePairWeights = rand(length(DNasePairFeats), 1);

for rr = 1:reclust;
    % Learn the regression for each module for each re-assignment of genes
    % to clusters
    
    fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
    fprintf('   Re-Clustering Iteration = %g/ %g\n', rr, reclust);
    fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
    
    nonEmptyModules = [];
    for k = 1:num_mod
        % Make a list of non-empty modules
        geneLocationsIndexes = find(ismember(assign, k));
        if isempty(geneLocationsIndexes)
            % The module is empty, so do not find its average
            continue
        end
        nonEmptyModules = vertcat(nonEmptyModules, k);
    end
    
    dist = Inf(num_genes,num_mod);
    modmean = Inf(num_mod,nsamples);
    
    for n = 1 : num_mod
        genesInMod = find(assign==n);
        if isempty(genesInMod)
            % There are no genes in the module, so do not find the module's
            % average
            continue
        end
        if length(genesInMod) > 1
            % There is more than 1 gene in the module, so find the average
            % of the genes in the modules
            modmean(n,:) = mean( expression(genesInMod,:),1 );
        else
            modmean(n,:) = expression(genesInMod,:);
        end
    end
    
    if isempty(DNaseRel)
        % No intial regulator relevances were given
        DNaseRel = zeros(nregs, num_genes, size(DNasePairFeats{1},3));
    end
    W = learn_weightsDNase (tf_expression, modmean, assign, DNaseRel, cellTypeMap); % Allowing only regulators to be features
    if rr == 1
        % Store the initial weights with no learned meta-priors
        params.W0 = W;
    end
    
    obj = calculate_objectiveDNase ( modmean, tf_expression, assign, nonEmptyModules, W, DNaseRel, cellTypeMap, DNasePairWeights, e );
    
    for iter = 1 : maxiters
        % Execute each iteration of Lirnet, and update the weights and
        % metapriors as necessary
        
        fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
        fprintf('   ITER = %g/ %g\n', iter, maxiters);
        fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
        
        DNasePairWeights = learn_regul_metapriorDNase(DNasePairFeats, cellTypeMap, DNasePairWeights,W,assign, nonEmptyModules, l1Upper,l1Lower,e);
        
        DNaseRel = calculate_regul_relevanceDNase (l1Upper, l1Lower, DNasePairFeats, DNasePairWeights); % Number of regulators by number of genes
        
        fprintf('learn weights...\n');
        W = learn_weightsDNase (tf_expression, modmean, assign, DNaseRel, cellTypeMap); % Allowing only regulators to be features
        
        preobj = obj;
        obj = calculate_objectiveDNase ( modmean, tf_expression, assign, nonEmptyModules, W, DNaseRel, cellTypeMap, DNasePairWeights, e );
     
        fprintf('\niter = %g\tobj = %g\tdel = %g\n\n',iter, obj, preobj-obj);
        
        if abs(preobj - obj) < 0.0001
            % The objective has not changed much, so stop
            break
        end
    end
    
    for ii = 1:num_mod
        % Iterate through modules and compute the regression error for each
        % module
        
        g = find(assign==ii);
        if isempty(g)
            % The module is empty, so do not consider it
            params.errors(rr,ii) = 0;
            continue
        end
        pred = tf_expression'*W(ii,:)'; % predictions
        params.errors(rr,ii) = sum ( sum ( (expression(g,:) -  repmat(pred',length(g),1)).^2)); %sse
        
        dist(:,ii) = sum ( (expression - repmat(pred',num_genes,1)).^2,2);
    
    end
    
    params.assign(:,rr) = assign;
    params.weights{rr} = W;
    params.DNasePairWeights{rr} = DNasePairWeights;
    
    if rr < reclust
        % Will be doing more re-clustering iterations, so move genes as
        % necessary
        % GENES ARE NOT ALLOWED TO FORM NEW MODULES
        [~,new_assign] = min(dist,[],2);
        nmoved = sum(new_assign~=assign);
        fprintf('genes moved %d\n',nmoved);

        assign = new_assign;
        
        if nmoved < 0.05*(num_genes);
            % Most genes are no longer switching clusters, so stop
            break
        end
    end
    
end

params.DNaseRel = DNaseRel;
