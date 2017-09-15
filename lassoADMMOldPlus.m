function [z, history] = lassoADMM(A, b, lambdaOne, rho, alpha)

% lasso  Solve lasso problem via ADMM
% Original code from S. Boyd
% Code modified by I. Kaplow
%
% Solves the following problem via ADMM:
% minimize || Ax - b||_2^2 + lambdaOne*|| x ||_1
% NOTE: NO 1/2
%
% The solution is returned in the vector x.
%
% history is a structure that contains the objective value, the primal and
% dual residual norms, and the tolerances for the primal and dual residual
% norms at each iteration.
%
% rho is the augmented Lagrangian parameter.
%
% alpha is the over-relaxation parameter (typical values for alpha are
% between 1.0 and 1.8).
%
% lambdaOne is the L1 penalty
%
% More information can be found in the paper linked at:
% http://www.stanford.edu/~boyd/papers/distr_opt_stat_learning_admm.html
%

t_start = tic;

%QUIET    = 0;
QUIET    = 1;
MAX_ITER = 1000;
ABSTOL   = 1e-6;
RELTOL   = 1e-4;

m = size(A, 1);
n = size(A, 2);

% save a matrix-vector multiply
Atb = A'*b;

x = zeros(n,1);
z = zeros(n,1);
u = zeros(n,1);

% cache the factorization
[L U] = factor(A, rho);

if ~QUIET
    fprintf('%3s\t%10s\t%10s\t%10s\t%10s\t%10s\n', 'iter', ...
      'r norm', 'eps pri', 's norm', 'eps dual', 'objective');
end

for k = 1:MAX_ITER

    % x-update
    q = 2*(Atb) + (2*rho)*(z - u);    % temporary value, ADDED 2* becuase no 1/2
    if( m >= n )    % if skinny
       x = U \ (L \ q);
    else            % if fat
       x = q/(2*rho) - (A'*(U \ ( L \ (A*q) )))/(2*rho)^2; % ADDED 2* because no 1/2
    end
    clear q

    % z-update with relaxation
    zold = z;
    x_hat = alpha*x + (1 - alpha)*zold;
    z = shrinkage(x_hat + (u/(2*rho)), lambdaOne/(2*rho)); % ADDED 2* because no 1/2, divided u by rho to match formula

    % u-update
    u = u + (rho*(x_hat - z)); % ADDED rho to match formula

    % diagnostics, reporting, termination checks
    history.objval(k)  = objective(A, b, lambdaOne, x, z);

    history.r_norm(k)  = norm(x - z);
    history.s_norm(k)  = norm(-rho*(z - zold));
    clear zold

    history.eps_pri(k) = sqrt(n)*ABSTOL + RELTOL*max(norm(x), norm(-z));
    history.eps_dual(k)= sqrt(n)*ABSTOL + RELTOL*norm(rho*u);

    if ~QUIET
        fprintf('%3d\t%10.4f\t%10.4f\t%10.4f\t%10.4f\t%10.2f\n', k, ...
            history.r_norm(k), history.eps_pri(k), ...
            history.s_norm(k), history.eps_dual(k), history.objval(k));
    end

    if (history.r_norm(k) < history.eps_pri(k) && ...
       history.s_norm(k) < history.eps_dual(k))
         break;
    end

end

if ~QUIET
    toc(t_start);
end

clear x
clear u

end

function p = objective(A, b, lambdaOne, x, z)
    p = ( sum((A*x - b).^2) + (lambdaOne.*abs(z)));
end

function z = shrinkage(x, kappa)
    z = max( 0, x - kappa ) - max( 0, -x - kappa );
end

function [L U] = factor(A, rho)
    [m, n] = size(A);
    if ( m >= n )    % if skinny
       L = chol( 2*(A'*A) + (2*rho)*speye(n), 'lower' ); % ADDED 2* because no 1/2
    else            % if fat
       L = chol( speye(m) + (1/rho)*((A*A')), 'lower' ); % Did not add 2* because no 1/2, but 2's for A and rho cancel
    end

    % force matlab to recognize the upper / lower triangular structure
    L = sparse(L);
    U = sparse(L');
end