function [Rtp,Rti,options] = linReachInner(sys,Rstart,params,options)
% linReachInner - computes an inner approximation of the reachable set
%    using the approach in [1], which is conceptually similar to [2]
%
% Syntax:
%    [Rtp,Rti,options] = linReachInner(sys,Rstart,params,options)
%
% Inputs:
%    sys - nonlinearSys object
%    Rstart - initial reachable set
%    params - model parameters
%    options - struct with algorithm settings
%
% Outputs:
%    Rtp - inner approximation of the time-point reachable set
%    Rti - inner approximation of the time-interval reachable set
%    options - struct with algorithm settings
%
% References: 
%   [1] M. Wetzlinger, A. Kulmburg, and M. Althoff. "Inner approximations
%       of reachable sets for nonlinear systems using the Minkowski
%       difference". IEEE Control Systems Letters, 2024.
%   [2] M. Althoff et al. "Reachability analysis of nonlinear systems with 
%       uncertain parameters using conservative linearization"
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contDynamics/linReach

% Authors:       Mark Wetzlinger
% Written:       17-December-2023
% Last update:   22-December-2023 (MW, major speed up)
%                28-November-2024 (MW, remove exponentialMatrix class)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% extract initial set and abstraction error
Rinit = Rstart.set;
abstrerr = Rstart.error;

% linearization -----------------------------------------------------------

% compute linearization point as center of interval enclosing the start set
% (replacing the Chebyshev center, which would result in a worse Lagrange
% remainder for very-non-box-like-shaped start sets; would also be slower)
Rinit_interval = interval(Rinit);
c_Rinit_interval = center(Rinit_interval);

% check whether inner approximation has become empty
if representsa_(Rinit_interval,'emptySet',1e-12)
    options.empty = true;
    Rtp = []; Rti = [];
    return
end

% formula: c_x + 1/2*dt*f(c_x,c_u)
options.linearizationPoint = c_Rinit_interval ...
    + 0.5*options.timeStep*sys.mFile(c_Rinit_interval,params.uTrans);

% linearize the nonlinear system
[sys,linsys,linOptions] = linearize(sys,Rinit,params,options);


% reachable set of linearized dynamics ------------------------------------
% time-point solution: linear map + constant input solution
% time-interval solution: compute outer approximation of the reachable set
% of the linearized system -> required to overestimate(!) Lagrange
% remainder over the time interval

% note: variables with suffix '__' are shifted by the linearization point

% translate Rinit by linearization point
Rinit__ = Rinit + (-sys.linError.p.x);

% constant input on right-hand side of differential inclusion
u = center(linOptions.uTrans);

% compute interval matrices F and G
[~,F,G] = taylorMatrices(linsys,options.timeStep,Inf);

% compute constant input solution (just a vector!)
Rcon = particularSolution_constant(linsys,u,options.timeStep,Inf);
% propagation matrix
eADeltat = getTaylor(linsys,'eAdt',struct('timeStep',options.timeStep));

% time-point solution (translated by linearization point)
Rlintp__ = eADeltat*Rinit__ + Rcon;

% box of translated initial set and curvature error
Rstart_interval__ = Rinit_interval + (-sys.linError.p.x);
C = F * zonotope(Rstart_interval__) + G*u;

% compute box enclosure of time-interval solution for abstraction error
% computation (uses box enclosure anyway) with formula
%    Rlinti = convHull_(Rstart,Rlintp) + C
% instead of interval(Rlintp), we use
%    eADeltat * Rstart_interval + Rcon
% because it is faster and the increase in the Lagrange remainder is small
Rlintp_outerapprox__ = eADeltat * Rstart_interval__ + Rcon;
Rlinti__ = convHull_(Rstart_interval__,Rlintp_outerapprox__ + interval(C),'exact');

% reachable set due to abstraction error ----------------------------------

% loop until the actual abstraction error is smaller than the estimated
% linearization error
perfIndCurr = Inf;

while perfIndCurr > 1
    % estimate the abstraction error
    appliedError = 1.1*abstrerr;
    Verror = zonotope(0*appliedError,diag(appliedError));
    RallError = interval(particularSolution_timeVarying(linsys,Verror,options.timeStep,Inf));

    % compute the time-interval reachable set including abstraction error
    % (represented as an interval since absterr_lin below converts to
    % interval anyway...)
    Rmax__ = Rlinti__ + RallError;
    % compute the abstraction error using the conservative
    % linearization approach described in [1]
    [trueError,VerrorDyn] = priv_abstrerr_lin(sys,zonotope(Rmax__),params,options);
    
    % compare linearization error with the maximum allowed error
    perfIndCurr = max(trueError./appliedError);
    abstrerr = trueError;
end

% compute the reachable set due to the linearization error
Rerror = particularSolution_timeVarying(linsys,VerrorDyn,options.timeStep,Inf);

% time-point and time-interval reachable sets -----------------------------

% time-point solution: Minkowski difference of reachable set for linearized
% dynamics (shifted by linearization point) and abstraction error
Rtp.set = minkDiff(Rlintp__ + sys.linError.p.x,Rerror);
% save abstraction error to use as initial guess in the next time step
Rtp.error = abstrerr;

% compute Minkowski difference of start set with curvature error and
% abstraction error for time-interval solution (only account for center of
% C since center(Rerror) is the origin)
Rstart_ti = minkDiff(minkDiff(Rinit,C) + 2*center(C),Rerror);
Rend_ti = minkDiff(Rtp.set,C) + 2*center(C);

% for the conversion to constrained zonotopes, any(!) enclosure of the
% polytope suffices, so we use Rmax shifted by the linearization point
Rmax = Rmax__ + sys.linError.p.x;
Rstart_ti_cZ = conZonotope(Rstart_ti,'exact:supportFunc',Rmax);
Rend_ti_cZ = conZonotope(Rend_ti,'exact:supportFunc',Rmax);
Rti = convHull_(Rstart_ti_cZ,Rend_ti_cZ,'exact:Raghuraman');

% ------------------------------ END OF CODE ------------------------------
