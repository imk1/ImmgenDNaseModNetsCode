function [old_dirs,old_stps,Hdiag] = lbfgsUpdateDNase(yMat,sMat,corrections,debug,old_dirs,old_stps,Hdiag)

y = reshape(yMat', size(yMat, 1) * size(yMat, 2), 1);
s = reshape(sMat', size(sMat, 1) * size(sMat, 2), 1);

if y'*s > 1e-10
    numCorrections = size(old_dirs,2);
    if numCorrections < corrections
        % Full Update
        old_dirs(:,numCorrections+1) = s;
        old_stps(:,numCorrections+1) = y;
    else
        % Limited-Memory Update
        old_dirs = [old_dirs(:,2:corrections) s];
        old_stps = [old_stps(:,2:corrections) y];
    end

    % Update scale of initial Hessian approximation
    Hdiag = (y'*s)/(y'*y);
else
    if debug
        fprintf('Skipping Update\n');
    end
end