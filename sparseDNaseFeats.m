function DNasePairFeats = sparseDNaseFeats(DNaseFeatFileNames)
% This function loads DNase feature matrices into Matlab and puts all of 
% the matrices into a cell array in which each matrix is sparse
% ASSUMES THAT ZEROING HAS ALREADY BEEN DONE

for i = 1:length(DNaseFeatFileNames)
    % Iterate through feature matrices and import all of them into Matlab
    i
    DNasePairFeats{i} = importdata(DNaseFeatFileNames{i});
    DNasePairFeats{i} = sparse(DNasePairFeats{i});
end

save /scr/ImmgenProject/MouseDNase/DNasePairFeatsSparse DNasePairFeats