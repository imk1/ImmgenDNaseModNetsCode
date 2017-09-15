function [DNasePairWeights,g,info] = spgDNase(funObj, funProj, DNasePairWeights, DNasePairFeats, W, assign, l1Upper, l1Lower, e, PWMScores, opts)
% SPG  Spectral Projected Gradient solver
%
%   [X,G,INFO] = SPG(FUNOBJ, FUNPROJ, XINIT, OPTS)
%
%   Solves the convex constrained problem
%
%        minimize  f(x)  subj to  x \in C
%
%   where C is an arbitrary convex set. The user-provided function
%
%        x = funProj(b)
%
%   computes the orthogonal projection of b onto C.
%
%   The following fields of the OPTS options structure are
%   recognized:
%
%   curvilinear	   0 or [1], disable or enable curvilinear backtracking
%   spectral       0 or [1], disable or enable spectral step length
%   maxiter        Maximum number of iterations
%   maxfunc        Maximum number of function evaluations
%   verbosity      0 = no output, [1] = all output
%   history        Number of previous function values of line search
%
% Apr 04 2007: First version, based on Birgin, Martinez, Raydan 2000
%              Michael P. Friedlander (mpf@cs.ubc.ca)
%
% $Id: spg.m 979 2008-06-01 09:00:56Z ewout78 $

if nargin < 11, opts = struct(); end;

if ~isfield(opts,'curvilinear'), opts.curvilinear = 0;               end; % Change back to 1!
if ~isfield(opts,'spectral'   ), opts.spectral = 1;               end;
if ~isfield(opts,'maxiter'    ), opts.maxiter     = 2;            end; % Change back to 1000!
if ~isfield(opts,'maxfunc'    ), opts.maxfunc     = 10*opts.maxiter; end;
if ~isfield(opts,'verbosity'  ), opts.verbosity   = 1;               end;
if ~isfield(opts,'history'    ), opts.history     = 10;              end;



% Parameters.
m       = opts.history;  % Number of previous function values.
optTol  = 1e-7;          % Optimality tol on projected gradient.
maxIts  = opts.maxiter;  % Maximum number of iterations.
maxFunc = opts.maxfunc;  % Maximum number of function evaluations.
stepMin = 1e-5;          % Minimum steplength.
stepMax = 1e+5;          % Maximum steplength.
fid     = 1;             % Output to screen.

spgFlag   = opts.spectral;
curvyFlag = opts.curvilinear;

logBody = ' %4i  %13.6e  %12.6e  %8.2e  %7.1e\n';
logHead = ' %4s  %13s  %12s  %8s  %7s\n';

if opts.verbosity
   fprintf(fid,logHead,'Iter','Objective','ProjGrad','gStep','lStep');
end

EXIT_OPTIMAL    = 0;
EXIT_ITERATIONS = 1;
EXIT_FUNCTIONS  = 2;

% Set history to one if spgFlag is not set
if spgFlag == 0
   m = 1;
end

% Initialize local variables.
tic;                   % Start your watches!
iter     = 0;
nObj     = 0;
nGrd     = 0;
lastFv   = -inf(m,1);   % Last m function values.
lStep    = 0;           % Linesearch step.
fHist    = zeros(1,maxIts+1);
normHist = zeros(1,maxIts+1);
evalHist = zeros(1,maxIts+1);

% Project the starting point and evaluate function and gradient.
DNasePairWeights = funProj(DNasePairWeights);
[f,g]            = funObj(DNasePairWeights, DNasePairFeats, W, assign, l1Upper, l1Lower, e, PWMScores);    nObj=nObj+1; nGrd=nGrd+1;
lastFv(1)        = f;
fBest            = f;
xBest            = DNasePairWeights;
fHist(1)         = f;
evalHist(1)      = nObj;


% Compute projected gradient direction and initial steplength.
d      = funProj(DNasePairWeights - g) - DNasePairWeights;
dNorm  = norm(d,inf);
gStep  = min( stepMax, max(stepMin, 1/dNorm) );
if spgFlag == 0, gStep = 1; end;

%----------------------------------------------------------------------
% MAIN LOOP.
%----------------------------------------------------------------------
while 1
    
    % Print to log.
    if opts.verbosity
       fprintf(fid,logBody,iter,f,dNorm,gStep,lStep);
    end
    
    % Test exit conditions.
    if dNorm < optTol
       stat = EXIT_OPTIMAL;
       break
    elseif iter >= maxIts
       stat = EXIT_ITERATIONS;
       break
    elseif nObj >= maxFunc
       stat = EXIT_FUNCTIONS;
       break
    end
       
    % Iterations begin here.
    iter = iter + 1;
    
    % Compute projected gradient direction.
    if ~curvyFlag
       d    = funProj(DNasePairWeights - gStep*g) - DNasePairWeights;
       gtd  = g'*d;
    end
    
    % Linesearch in the direction d.
    if ~curvyFlag
       [fNew,xNew,lStep,lnObj] = spgLine(f,DNasePairWeights,d,gtd,max(lastFv),funObj, DNasePairFeats, W, assign, l1Upper, l1Lower, e, PWMScores);
    else
       [fNew,xNew,lStep,lnObj] = spgLineCurvy(DNasePairWeights,gStep*g,max(lastFv),funObj,funProj, DNasePairFeats, W, assign, l1Upper, l1Lower, e, PWMScores);
    end
 
    lastFv(mod(iter,m)+1) = fNew;
    nObj = nObj + lnObj;
    
    % Update to best iterate.
    if fBest > fNew
       fBest = fNew;
       xBest = xNew;
    end
   
    % Update gradient and compute new step length.
    [fNew,gNew] = funObj(xNew, DNasePairFeats, W, assign, l1Upper, l1Lower, e, PWMScores); nObj=nObj+1; nGrd=nGrd+1;
    [dummy,groupNorm] = funProj(xNew);
    fHist(iter+1)     = fNew;
    normHist(iter+1)  = groupNorm;
    evalHist(iter+1)  = nObj;

    s = xNew - DNasePairWeights;
    y = gNew - g;
    DNasePairWeights = xNew;
    g = gNew;
    f = fNew;

    % Quantities needed to test optimality.
    d = funProj(DNasePairWeights - g) - DNasePairWeights;
    dNorm = norm(d,inf);

    sts = s'*s;
    sty = s'*y;
    if sty  <= 0
       gStep  = stepMax;
    else
       gStep  = min( stepMax, max(stepMin, sts/sty) );
    end
    
    if spgFlag == 0, gStep = 1; end;
    
end % while 1

% Restore best solution (if needed).
if f > fBest
   if opts.verbosity
      fprintf(fid,' Restoring best iterate.\n');
   end
   DNasePairWeights = xBest;
   [f,g] = funObj(DNasePairWeights, DNasePairFeats, W, assign, l1Upper, l1Lower, e, PWMScores); nObj=nObj+1; nGrd=nGrd+1;
end

info.stat     = stat;
info.nObj     = nObj;
info.nGrd     = nGrd;
info.iter     = iter;
info.fHist    = fHist(1:iter+1);
info.normHist = normHist(1:iter+1);
info.evalHist = evalHist(1:iter+1);

switch(stat)
   case EXIT_OPTIMAL
      info.statstr = 'Optimal solution found';

   case EXIT_ITERATIONS
      info.statstr = 'Too many iterations';

   case EXIT_FUNCTIONS
      info.statstr = 'Too many function evaluations';

   otherwise
      error('Unknown termination condition\n');
end

if opts.verbosity
   fprintf(fid,'\n');
   fprintf(fid,' EXIT -- %s\n',info.statstr)
   fprintf(fid,'\n');
   fprintf(fid,' Function evaluations: %6i  '  ,nObj);
   fprintf(fid,' Total time (sec):     %6.1f\n',toc );
   fprintf(fid,' Gradient evaluations: %6i\n'  ,nGrd);
end

end % function spg

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fNew,xNew,step,nObj,info] = spgLine(f,x,d,gtd,fMax,funObj, DNasePairFeats, W, assign, l1Upper, l1Lower, e, PWMScores)
% Nonmonotone linesearch.

EXIT_CONVERGED  = 1;
EXIT_ITERATIONS = 2;
maxIts = 10;
step   = 1;
iter   = 0;
nObj   = 0;
gamma  = 1e-4;

while 1

    % Evaluate trial point and function value.
    xNew = x + step*d;
    fNew = funObj(xNew, DNasePairFeats, W, assign, l1Upper, l1Lower, e, PWMScores);
    nObj = nObj + 1;

    % Check exit conditions.
    if fNew < fMax + gamma*step*gtd  % Sufficient descent condition.
       info = EXIT_CONVERGED;
       break
    elseif  iter >= maxIts           % Too many linesearch iterations.
       info = EXIT_ITERATIONS;
       break
    end
    
    % New linesearch iteration.
    iter = iter + 1;
    
    % Safeguarded quadratic interpolation.
    if step <= 0.1
       step  = step / 2;
    else
       tmp = (-gtd*step^2) / (2*(fNew-f-step*gtd));
       if tmp < 0.1 || tmp > 0.9*step || isnan(tmp)
          tmp = step / 2;
       end
       step = tmp;
    end
    
end % while 1

end % function spgLine


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fNew,xNew,step,nObj,info] = ...
    spgLineCurvy(x,g,fMax,funObj,funProj, DNasePairFeats, W, assign, l1Upper, l1Lower, e, PWMScores)
% Projected backtracking linesearch.
% On entry,
% g  is the (possibly scaled) steepest descent direction.

EXIT_CONVERGED  = 0;
EXIT_ITERATIONS = 1;
EXIT_NODESCENT  = 2;
gamma  = 1e-4;
maxIts = 10;
step   =  1;
sNorm  =  0;
scale  =  1;      % Safeguard scaling.  (See below.)
nSafe  =  0;      % No. of safeguarding steps.
nObj   =  0;
iter   =  0;
debug  =  false;  % Set to true to enable log.
n      =  length(x);

if debug
   fprintf(' %5s  %13s  %13s  %13s  %8s\n',...
           'LSits','fNew','step','gts','scale');  
end
   
while 1

    % Evaluate trial point and function value.
    xNew     = funProj(x - step*scale*g);
    fNew     = funObj(xNew, DNasePairFeats, W, assign, l1Upper, l1Lower, e, PWMScores);
    nObj     = nObj + 1;
    s        = xNew - x;
    gts      = scale * g' * s;
    if gts >= 0 % Should we check real and complex parts individually?
       info = EXIT_NODESCENT;
       break
    end

    if debug
       fprintf(' LS %2i  %13.7e  %13.7e  %13.6e  %8.1e\n',...
               iter,fNew,step,gts,scale);
    end
    
    % 03 Aug 07: If gts is complex, then should be looking at -abs(gts).
    if fNew < fMax - gamma*step*abs(gts)  % Sufficient descent condition.
       info = EXIT_CONVERGED;
       break
    elseif iter >= maxIts                 % Too many linesearch iterations.
       info = EXIT_ITERATIONS;
       break
    end
    
    % New linesearch iteration.
    iter = iter + 1;
    step = step / 2;

    % Safeguard: If stepMax is huge, then even damped search
    % directions can give exactly the same point after projection.  If
    % we observe this in adjacent iterations, we drastically damp the
    % next search direction.
    % 31 May 07: Damp consecutive safeguarding steps.
    sNormOld  = sNorm;
    sNorm     = norm(s) / sqrt(n);
    %   if sNorm >= sNormOld
    if abs(sNorm - sNormOld) <= 1e-6 * sNorm
       gNorm = norm(g) / sqrt(n);
       scale = sNorm / gNorm / (2^nSafe);
       nSafe = nSafe + 1;
    end
    
end % while 1

end % function spgLineCurvy
