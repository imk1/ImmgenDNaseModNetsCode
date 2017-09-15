function obj = calculate_objectiveDNase ( modmean, expression, W, DNAseRel, d, e )

% Code written by Su-In Lee
% Modified by Irene Kaplow

% Calculates the full objective for DNase Lirnet

ll = sum( sum( (modmean - W * expression).^2 ) );
obj = ll + sum( sum( abs(W) .* DNAseRel' ) ) + (d*sum(sum(W.^2))) + (e*sum(DNAseRel.^2));