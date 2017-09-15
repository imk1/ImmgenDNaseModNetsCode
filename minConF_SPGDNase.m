function [x,f,funEvals,projects] = minConF_SPGDNase(funObj,x,funProj, options)
% function [x,f] = minConF_SPG(funObj,x,funProj,options)
%
% Function for using Spectral Projected Gradient to solve problems of the form
%   min funObj(x) s.t. x in C
%
%   @funObj(x): function to minimize (returns gradient as second argument)
%   @funProj(x): function that returns projection of x onto C
%
%   options:
%       verbose: level of verbosity (0: no output, 1: final, 2: iter (default), 3:
%       debug)
%       optTol: tolerance used to check for progress (default: 1e-6)
%       maxIter: maximum number of calls to funObj (default: 500)
%       numDiff: compute derivatives numerically (0: use user-supplied
%       derivatives (default), 1: use finite differences, 2: use complex
%       differentials)
%       suffDec: sufficient decrease parameter in Armijo condition (default
%       : 1e-4)
%       interp: type of interpolation (0: step-size halving, 1: quadratic,
%       2: cubic)
%       memory: number of steps to look back in non-monotone Armijo
%       condition
%       useSpectral: use spectral scaling of gradient direction (default:
%       1)
%       curvilinear: backtrack along projection Arc (default: 0)
%       testOpt: test optimality condition (default: 1)
%       feasibleInit: if 1, then the initial point is assumed to be
%       feasible
%       bbType: type of Barzilai Borwein step (default: 1)
%
%   Notes: 
%       - if the projection is expensive to compute, you can reduce the
%           number of projections by setting testOpt to 0

% Set Parameters
if nargin < 11
    options = [];
end
[verbose,optTol,maxIter,suffDec,interp,memory,useSpectral,curvilinear,feasibleInit,testOpt,bbType] = ...
    myProcessOptions(...
    options,'verbose',2,'optTol',1e-6,'maxIter',500,'suffDec',1e-4,...
    'interp',2,'memory',10,'useSpectral',1,'curvilinear',0,'feasibleInit',0,...
    'testOpt',1,'bbType',1);
% Assumes that numDiff is 0

% Output Log
if verbose >= 2
    if testOpt
        fprintf('%10s %10s %10s %15s %15s %15s\n','Iteration','FunEvals','Projections','Step Length','Function Val','Opt Cond');
    else
        fprintf('%10s %10s %10s %15s %15s\n','Iteration','FunEvals','Projections','Step Length','Function Val');
    end
end

xCol = reshape(x', size(x, 1) * size(x, 2), 1);

% Make objective function (if using numerical derivatives)
funEvalMultiplier = 1;

% Evaluate Initial Point
if ~feasibleInit
    xCol = funProj(xCol);
end
[f,g] = funObj(x);
projects = 1;
funEvals = 1;

% Optionally check optimality
if testOpt
    projects = projects+1;
    if sum(abs(funProj(xCol-g)-xCol), 2) < optTol
        if verbose >= 1
        fprintf('First-Order Optimality Conditions Below optTol at Initial Point\n');
        end
        return;
    end
end

i = 1;
while funEvals <= maxIter

    % Compute Step Direction
    if i == 1 || ~useSpectral
        alpha = 1;
    else
        y = g-g_old;
        s = xCol-xOldCol;
        if bbType == 1
            alpha = (s'*s)/(s'*y);
        else
            alpha = (s'*y)/(y'*y);
        end
        if alpha <= 1e-10 || alpha > 1e10
            alpha = 1;
        end
    end
    d = -alpha*g;
    f_old = f;
    xOldCol = xCol;
    g_old = g;
    
    % Compute Projected Step
    if ~curvilinear
        d = funProj(xCol+d)-xCol;
        projects = projects+1;
    end

    % Check that Progress can be made along the direction
    gtd = g'*d;
    if gtd > -optTol
        if verbose >= 1
            fprintf('Directional Derivative below optTol\n');
        end
        break;
    end

    % Select Initial Guess to step length
    if i == 1
        t = min(1,1/sum(abs(g)));
    else
        t = 1;
    end

    % Compute reference function for non-monotone condition

    if memory == 1
        funRef = f;
    else
        if i == 1
            old_fvals = repmat(-inf,[memory 1]);
        end

        if i <= memory
            old_fvals(i) = f;
        else
            old_fvals = [old_fvals(2:end);f];
        end
        funRef = max(old_fvals);
    end

    % Evaluate the Objective and Gradient at the Initial Step
    if curvilinear
        xNewCol = funProj(xCol + t*d);
        projects = projects+1;
    else
        xNewCol = xCol + t*d;
    end
    x_new = reshape(xNewCol, size(x,2), size(x,1))';
    [f_new,g_new] = funObj(x_new);
    funEvals = funEvals+1;

    % Backtracking Line Search
    lineSearchIters = 1;
    maxLineSearchIters = 100;
    while (f_new > funRef + suffDec*g'*(xNewCol-xCol) || ~isLegal(f_new)) && (lineSearchIters <= maxLineSearchIters)
        temp = t;
        if interp == 0 || ~isLegal(f_new)
            if verbose == 3
                fprintf('Halving Step Size\n');
            end
            t = t/2;
        elseif interp == 2 && isLegal(g_new)
            if verbose == 3
                fprintf('Cubic Backtracking\n');
            end
            t = polyinterp([0 f gtd; t f_new g_new'*d]);
        elseif lineSearchIters < 2 || ~isLegal(f_prev)
            if verbose == 3
                fprintf('Quadratic Backtracking\n');
            end
            t = polyinterp([0 f gtd; t f_new sqrt(-1)]);
        else
            if verbose == 3
                fprintf('Cubic Backtracking on Function Values\n');
            end
            t = polyinterp([0 f gtd; t f_new sqrt(-1);t_prev f_prev sqrt(-1)]);
        end

        % Adjust if change is too small
        if t < temp*1e-3
            if verbose == 3
                fprintf('Interpolated value too small, Adjusting\n');
            end
            t = temp*1e-3;
        elseif t > temp*0.6
            if verbose == 3
                fprintf('Interpolated value too large, Adjusting\n');
            end
            t = temp*0.6;
        end

        % Check whether step has become too small
        if sum(abs(t*d)) < optTol || t == 0
            if verbose == 3
                fprintf('Line Search failed\n');
            end
            t = 0;
            f_new = f;
            g_new = g;
            break;
        end

        % Evaluate New Point
        f_prev = f_new;
        t_prev = temp;
        if curvilinear
            xNewCol = funProj(xCol + t*d);
            projects = projects+1;
        else
            xNewCol = xCol + t*d;
        end
        x_new = reshape(xNewCol, size(x,2), size(x,1))';
        [f_new,g_new] = funObj(x_new);
        funEvals = funEvals+1;
        lineSearchIters = lineSearchIters+1;

    end

    % Take Step
    x = x_new;
    xCol = xNewCol;
    f = f_new;
    g = g_new;

    if testOpt
        optCond = sum(abs(funProj(xCol-g)-xCol));
        projects = projects+1;
    end

    % Output Log
    if verbose >= 2
        if testOpt
            fprintf('%10d %10d %10d %15.5e %15.5e %15.5e\n',i,funEvals*funEvalMultiplier,projects,t,f,optCond);
        else
            fprintf('%10d %10d %10d %15.5e %15.5e\n',i,funEvals*funEvalMultiplier,projects,t,f);
        end
    end

    % Check optimality
    if testOpt
        if optCond < optTol
            if verbose >= 1
            fprintf('First-Order Optimality Conditions Below optTol\n');
            end
            break;
        end
    end

    if sum(abs(t*d)) < optTol
        if verbose >= 1
            fprintf('Step size below optTol\n');
        end
        break;
    end

    if abs(f-f_old) < optTol
        if verbose >= 1
            fprintf('Function value changing by less than optTol\n');
        end
        break;
    end

    if funEvals*funEvalMultiplier > maxIter
        if verbose >= 1
            fprintf('Function Evaluations exceeds maxIter\n');
        end
        break;
    end

    i = i + 1;
end
end