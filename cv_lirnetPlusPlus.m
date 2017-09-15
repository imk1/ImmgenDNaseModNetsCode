function [error,rsquare,rsquareModule,weights,assign] = cv_lirnetPlusPlus(X,Y,l1,l2,nfolds)

% Code written by Sara Mostafavi
% Modified by Irene Kaplow

% This code learns module networks using elastic net regularization and
% nfolds-fold cross-validation

% Added -- Assumes that X and Y are genes x cell types
X = X';
Y = Y';

[n,d] = size(X);
rand('seed',0);

pp = randperm(n);

init = 0;

nn = round(n/nfolds);

for ii = 1:nfolds
    % Make each validation set
    idx{ii} = pp(init+1:min(n,init+nn));
    init = init+nn;
end

for ii = 1:nfolds
    % Do training and cross-validation with each validation set
    
    diff = setdiff(1:n, idx{ii});
    %diff = setdiff(1:n,ii);
    
    Xtr = X(diff,:);
    Ytr = Y(diff,:);    
    Xtrn = whiten(Xtr')';
    Ytrn = whiten(Ytr')';
    
    Xtest = X(idx{ii}, :);
    Ytest = Y(idx{ii}, :);
    %Xtest = X(ii,:);
    %Ytest = Y(ii,:);
    
    Xtest = (Xtest-ones(size(Xtest, 1),1)*mean(Xtr))./(ones(size(Xtest, 1),1)*sqrt(sum(Xtr.^2)));
    Ytest = (Ytest-ones(size(Ytest, 1),1)*mean(Ytr))./(ones(size(Ytest, 1),1)*sqrt(sum(Ytr.^2)));
    %Xtest = (Xtest-mean(Xtr))./sqrt(sum(Xtr.^2));
    %Ytest = (Ytest-mean(Ytr))./sqrt(sum(Ytr.^2));
    
    ii
    for jj = 1:length(l1)
        % Train with each (l1, l2) parameter pair
        
        params = lirnet_flat_fullPlusPlus(Xtrn',Ytrn',l1(jj),l2(jj)); % X is all expressions, Y is TF expressions
        pred = Ytest*params.weights{end};
        %params = lirnet_flat_fullPlus(Ytrn',Xtrn',num_mod,l1(jj),l2(jj));
        %pred = Xtest*params.weights{end};
        
        pred = pred(:,params.assign(:,end));
        
        error{ii,jj} = ((pred-Xtest).^2)./(Xtest.^2)
        rsquare(ii,jj) = (corr(pred(:),Xtest(:))).^2;
        %error{ii,jj} = ((pred-Ytest).^2)./(Ytest.^2)
        %rsquare(ii,jj) = (corr(pred(:),Ytest(:))).^2;
        
        indexes = unique(params.assign(:,end));
        idxPlus = zeros(length(params.assign(:,end)), 1);
        for i = 1:length(params.assign(:,end))
            index = find(indexes == params.assign(i,end));
            idxPlus(i) = index;
        end
        for k = 1:max(idxPlus)
            % Iterate through modules and find the r squared for each
            % module
            modLocations = find(idxPlus == k);
            if length(modLocations) > 0
                % There is at least 1 gene in the current module
                predModule = pred(:,modLocations);
                XtestModule = Xtest(:,modLocations);
                rsquareModule(ii,jj,k) = (corr(predModule(:),XtestModule(:))).^2;
            else
                rsquareModule(ii,jj,k) = NaN;
            end
        end
        
        weights{ii,jj} = params.weights{end};
        
        assign{ii,jj} = params.assign(:,end);
    end
end
    