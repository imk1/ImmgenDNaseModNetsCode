function DNasePairWeights = learn_regul_metapriorDNase (DNasePairFeats, DNasePairWeights, W, assign, l1Upper, l1Lower, e, PWMScores)

% Code written by Su-In Lee
% Modified by Irene Kaplow

% Calculates the meta-priors for the DNase features

fprintf('BOF::learn_regul_metaprior\n');

u = abs(W)';

alpha = 0.3;
beta = 0.8;
mu = 10;
T = 1; % Maybe look into duality gap and change

ngenes = size( DNasePairFeats{1},2 );

for barrierIter = 1:100
    % Iterate through barrier method iterations until convergence or until
    % the maximum number of iterations has been done
    for iter = 1 : 100
        % Iterate through gradient descent iterations until convergence or until
        % the maximum number of iterations has been done
        % Use barrier method to constrain feature weights to be < 0 [add
        % -log(-x) to objective, multiply original objective by T > 0, and
        % multiply T by mu > 1 on every iteration

        DNaseRel = calculate_regul_relevanceDNase (l1Upper, l1Lower, DNasePairFeats, DNasePairWeights, PWMScores);
        DNaseRelAverage = zeros(size(u));
        for k = 1:size(u, 2)
            % Compute the average relevance for each regulator across modules
            geneLocationsIndexes = find(ismember(assign, k));
            DNaseRelAverage(:,k) = mean(DNaseRel(:,geneLocationsIndexes), 2)'; % Using 1 l1 penalty that is the average of all of the l1 penalties for the genes in the module
        end
        preobj = e * (DNasePairWeights')*DNasePairWeights + sum( sum( u .* DNaseRelAverage ) );

        pgrad = 2 * e * DNasePairWeights;

        for m = 1 : ngenes;
            % Compute the gradient with respect to the weights for each gene
            rel = zeros( size(DNasePairFeats{1},1), 1);
            for k = 1 : length(DNasePairFeats)
                % Iterate through pairwise features and add them appropriately
                % Assumes that pairwise features are 0 for DNases that are not
                % relevant to the gene
                rel = rel + DNasePairFeats{k}(:,m) * DNasePairWeights(k);
            end
            relTF = PWMScores * rel;
            z = computeSigmoid( relTF );
            for k = 1 : length(DNasePairFeats)
                % Compute the gradient with respect to the weight for each feature
                preSumFeats = PWMScores;
                for i = 1:size(PWMScores, 1)
                    % Multiply each PWMScore by the corresponding DNaseFeature
                    preSumFeats(i,:) = preSumFeats(i,:) .* DNasePairFeats{k}(:,m)';
                end
                sumFeats = sum(preSumFeats, 2);
                pgrad(k) = pgrad(k) + ((l1Lower-l1Upper) * ( ( z .* (1-z) .* u(:,assign(m)) )'*sumFeats));
                % Account for the Laplacian normalization constant
                pgrad(k) = pgrad(k) - (((l1Lower-l1Upper) * ( ( z .* (1-z) .* u(:,assign(m)) ).*sumFeats))' * (ones(size(z)) ./ ((l1Lower - l1Upper) * z + l1Upper)));
                % Multiply by T and add the derivative of 
                % -log(-feature weight)
                pgrad(k) = (T*pgrad(k)) - (1/(DNasePairWeights(k)));
            end
        end

        grad2 = pgrad'*pgrad;

        t = 1;

        for biter = 1 : 100
            % Iterate through backtracking line search iterations to choose the step size

            new_pb = DNasePairWeights - t * pgrad;

            new_a = calculate_regul_relevanceDNase ( l1Upper, l1Lower, DNasePairFeats, new_pb, PWMScores);
            new_aAverage = zeros(size(u));
            for k = 1:size(u, 2)
                % Compute the average relevance for each regulator across modules
                geneLocationsIndexes = find(ismember(assign, k));
                new_aAverage(:,k) = mean(new_a(:,geneLocationsIndexes), 2)'; % Using 1 l1 penalty that is the average of all of the l1 penalties for the genes in the module
            end
            postobj = e * (new_pb')*new_pb + sum( sum( u .* new_aAverage ) );

            if( postobj < preobj - alpha * t * grad2 )
                % Step-size leads to a sufficient decrease, so stop
                break;
            else
                t = t * beta;
            end
        end

        if( postobj<preobj )
            % Step has led to decrease, so take the step
            DNasePairWeights = new_pb;
        end

        fprintf('iter = %g\t%g\t%g\tdel = %g\n',iter, preobj, postobj, preobj-postobj);

        if( abs( preobj-postobj )<0.001 )
            % Objective is no longer changing much, so stop
            break;
        end
    end
    if length(DNasePairWeights)/T < 0.0001
        % T is no longer changing much, so stop
        break
    end
    T = mu * T;
end

fprintf('EOF::learn_marker_metaprior\n\n');
