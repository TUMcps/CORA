function timeInt = priv_reach_backward_AE_timeinterval(linsys,params,options)
% priv_reach_backward_AE_timeinterval - computes an outer approximation of the
%    time-interval minimal backward reachable set defined as
%       R(-tau) := { x_0 in R^n | exists u() forall w() exists t in tau:
%                                 x(t; x_0, u(), w()) in Xend }
%    via [1, Sec. V-B]
%
% Syntax:
%    timeInt = priv_reach_backward_AE_timeinterval(linsys,params,options)
%
% Inputs:
%    linsys - linearSys object
%    params - model parameters
%    options - reachability settings
%
% Outputs:
%    timeInt - struct of time-point solution
%
% References:
%    [1] M. Wetzlinger and M. Althoff, "Backward Reachability Analysis with
%        State Constraints for Linear Systems using Set Propagation", 2023.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%-
% See also: reachBackward

% Authors:       Mark Wetzlinger
% Written:       12-July-2023
% Last update:   28-November-2024 (MW, remove exponentialMatrix class)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% shorter variable names...
t0 = params.tStart;
dt = options.timeStep;
% time period and number of steps
tVec = t0:dt:params.tFinal;
steps = length(tVec) - 1;

% dimension
n = linsys.nrOfDims;

% init time-inverted dynamics (note: system is put in canonical form in
% reachBackward, so B = E = identity)
linsys_backward = linearSys(-linsys.A);
linsys_backward.taylor = taylorLinSys(linsys_backward.A);

% initialize struct for backward reachable set
% array of time intervals
timeIntArray = cellfun(@(x) interval(x(1),x(2)), ...
    num2cell([-tVec(2:end)',-tVec(1:end-1)'],2),'UniformOutput',false);
timeInt.set = cell(steps,1);
timeInt.time = timeIntArray;

% convert target set to a constrained zonotope
Rend_cZ = conZonotope(params.R0);
[Rend_isInterval,Rend_I] = representsa_(params.R0,'interval',1e-10);

% compute correction matrix for that state for given time step size
F = priv_correctionMatrixState(linsys_backward,dt,Inf);
if Rend_isInterval
    curvatureError_0_dt = F * zonotope(Rend_I);
    c_curvatureError_0_dt = center(curvatureError_0_dt);
    G_curvatureError_0_dt = generators(curvatureError_0_dt);
end

% init direction for support function evaluation (we need back-propagation
% in both forward (e^At) and backward (e^-At) perspective, for particular
% solutions and target set, respectively)
Norig = [eye(n), -eye(n)];
numDir = size(Norig,2);

% pre-compute particular solutions until start of considered time interval
% ...use 'closest' time step size to the one provided that partitions the
%    interval [0,params.tStart] into a whole number of steps
if withinTol(t0,0)
    steps_start = 0;
    ZW_outer_0_t = zeros(n,1);
    ZW_outer_0_t_sF = zeros(numDir,1+steps);
    c_ZW_outer_0_t_sF = zeros(numDir,1+steps);
    ZU_inner_0_t_sF = zeros(numDir,1+steps);

else
    % integer number of steps from 0 to tStart (search for a similar time
    % step size to the given time step size)
    steps_start = ceil(t0/dt);
    dt_ = t0/steps_start;

    % compute inner and outer approximations of particular solution [0,dt_]
    ZW_outer_0_dt = particularSolution_timeVarying(linsys,-params.W,dt_,Inf);
    c_ZW_outer_0_dt = center(ZW_outer_0_dt);
    G_ZW_outer_0_dt = generators(ZW_outer_0_dt);
    ZU_inner_0_dt = particularSolution_constant(linsys,-params.U,dt_,Inf);
    c_ZU_inner_0_dt = center(ZU_inner_0_dt);
    G_ZU_inner_0_dt = generators(ZU_inner_0_dt);

    % init particular solutions
    ZW_outer_0_t = zeros(n,1);
    ZW_outer_0_t_sF = zeros(numDir,1+steps_start+steps);
    c_ZW_outer_0_t_sF = zeros(numDir,1+steps_start+steps);
    ZU_inner_0_t_sF = zeros(numDir,1+steps_start+steps);

    eAdt = expm(linsys_backward.A*dt_);
    eAtk = eAdt;

    % propagate particular solutions until start time
    for k=1:steps_start
        % propagate explicit particular solutions
        ZW_outer_0_t = ZW_outer_0_t + eAtk * ZW_outer_0_dt;

        % propagate direction
        Nminus = eAtk * Norig;

        % loop over all directions
        ZW_outer_0_t_sF(:,k+1) = ZW_outer_0_t_sF(:,k) ...
            + Nminus' * c_ZW_outer_0_dt + sum(abs(Nminus' * G_ZW_outer_0_dt),2);
        c_ZW_outer_0_t_sF(:,k+1) = c_ZW_outer_0_t_sF(:,k) ...
            + Nminus' * c_ZW_outer_0_dt;
        ZU_inner_0_t_sF(:,k+1) = ZU_inner_0_t_sF(:,k) ...
            - Nminus' * c_ZU_inner_0_dt + sum(abs(-Nminus' * G_ZU_inner_0_dt),2);

        % propagate propagation matrix
        eAtk = eAtk * eAdt;
    end
end

% compute inner- and outer-approximations of particular solution [0,dt]
% (we may have a different time step size to before)
if steps_start == 0 || ~withinTol(dt_,dt)
    ZW_outer_0_dt = particularSolution_timeVarying(linsys,-params.W,dt,Inf);
    if isa(ZW_outer_0_dt,'contSet')
        c_ZW_outer_0_dt = center(ZW_outer_0_dt);
        G_ZW_outer_0_dt = generators(ZW_outer_0_dt);
    elseif isnumeric(ZW_outer_0_dt)
        c_ZW_outer_0_dt = ZW_outer_0_dt;
        G_ZW_outer_0_dt = zeros(n,1);
    end
    ZU_inner_0_dt = particularSolution_constant(linsys,-params.U,dt,Inf);
    if isa(ZU_inner_0_dt,'contSet')
        c_ZU_inner_0_dt = center(ZU_inner_0_dt);
        G_ZU_inner_0_dt = generators(ZU_inner_0_dt);
    elseif isnumeric(ZU_inner_0_dt)
        c_ZU_inner_0_dt = ZU_inner_0_dt;
        G_ZU_inner_0_dt = zeros(n,1);
    end
end


% compute explicit outer-approximation of maximal backward reachable set
% (no disturbances yet)
BRSE = cell(steps,1);

% pre-compute propagation matrices
eAtk = computeField(linsys_backward.taylor,'eAdt',struct('timeStep',t0));
if steps_start == 0 || ~withinTol(dt,dt_)
    eAdt = expm(linsys_backward.A*dt);
end

% support function evaluation of target set at start time
Xend_sF = zeros(numDir,1+steps);
Nminus = readFieldForTimeStep(linsys_backward.taylor,'eAdt',t0)' * Norig;
for j=1:numDir
    if Rend_isInterval
        Xend_sF(j,1) = supportFunc_(Rend_I,Nminus(:,j),'upper');
    else
        Xend_sF(j,1) = supportFunc_(Rend_cZ,Nminus(:,j),'upper');
    end
end

% offset for halfspaces
poly_offset = [Xend_sF(:,1), zeros(numDir,steps)];
BRSE_sF_lb = zeros(numDir,steps);

% loop over all steps from tStart to tFinal
for k=1:steps

    % propagate propagation matrix
    eAtk = eAtk * eAdt;
    insertFieldTimeStep(linsys_backward.taylor,'eAdt',eAtk,t0+dt*k);

    % propagate particular solutions
    ZW_outer_0_t = ZW_outer_0_t + ...
        readFieldForTimeStep(linsys_backward.taylor,'eAdt',t0+(k-1)*dt) * ZW_outer_0_dt;

    % propagate support function of particular solutions
    ZW_outer_0_t_sF(:,k+1) = ZW_outer_0_t_sF(:,k) ...
        + Nminus' * c_ZW_outer_0_dt + sum(abs(Nminus' * G_ZW_outer_0_dt),2);
    c_ZW_outer_0_t_sF(:,k+1) = c_ZW_outer_0_t_sF(:,k) ...
        + Nminus' * c_ZW_outer_0_dt;
    ZU_inner_0_t_sF(:,k+1) = ZU_inner_0_t_sF(:,k) ...
        - Nminus' * c_ZU_inner_0_dt + sum(abs(-Nminus' * G_ZU_inner_0_dt),2);

    % back-propagate directions: e^(At_k), e^(-At_k)
    Nminus = readFieldForTimeStep(linsys_backward.taylor,'eAdt',t0+k*dt)' * Norig;

    % loop over all directions
    for j=1:numDir
        % propagate support function of back-propagated target set
        if Rend_isInterval
            Xend_sF(j,k+1) = supportFunc_(Rend_I,Nminus(:,j),'upper');
        else
            Xend_sF(j,k+1) = supportFunc_(Rend_cZ,Nminus(:,j),'upper');
        end
    end

    % compute intermediate backward reachable set R_E(-tau_k)

    % explicit and implicit computation of curvature error
    if Rend_isInterval
        % faster version
        curvatureError = readFieldForTimeStep(linsys_backward.taylor,'eAdt',t0+(k-1)*dt) * curvatureError_0_dt;
        curvatureError_sF = Nminus' * c_curvatureError_0_dt ...
            - sum(abs(Nminus' * G_curvatureError_0_dt),2);

    else
        curvatureError = F * ...
            (readFieldForTimeStep(linsys_backward.taylor,'eAdt',t0+(k-1)*dt)*Rend_cZ);
        curvatureError_sF = zeros(numDir,1);
        for j=1:numDir
            curvatureError_sF(j) = ...
                -supportFunc_(curvatureError,-Norig(:,j),'upper');
        end
    end

    % update offset of intersecting polytope
    poly_offset(:,k+1) = max(Xend_sF(:,k:k+1),[],2) ...
        + ZW_outer_0_t_sF(:,k+1) - ZU_inner_0_t_sF(:,k);

    % explicit set for all-zero input trajectory
    BRSE{k} = convHull(readFieldForTimeStep(linsys_backward.taylor,'eAdt',t0+(k-1)*dt)*Rend_cZ,...
        readFieldForTimeStep(linsys_backward.taylor,'eAdt',t0+k*dt)*Rend_cZ) + curvatureError + ZW_outer_0_t;
    
    % lower bound of support function in each direction
    BRSE_sF_lb(:,k) = ...
        -max([Xend_sF(n+1:end,k:k+1); Xend_sF(1:n,k:k+1)],[],2) ...
        + curvatureError_sF ...
        + [ZW_outer_0_t_sF(n+1:end,k+1); ZW_outer_0_t_sF(1:n,k+1)] ...
        - 2*c_ZW_outer_0_t_sF(:,k+1);

end

timeInt.set = cell(steps,1);
% intersect backward reachable set for zero-input-trajectory with support
% function of inner approximations for other input trajectories
for k=1:steps
    timeInt.set{k,1} = aux_intersection(BRSE{k}, Norig, ...
        poly_offset(:,k+1), BRSE_sF_lb(:,k) - 1);
end

end


% Auxiliary functions -----------------------------------------------------

function cZ_and = aux_intersection(cZ,C,d_ub,d_lb)
% evaluate intersection of forward reachable set of time-inverted dynamics
% with center input trajectory aka BRSE (represented as a conZonotope) and
% same set using different extreme input trajectories (represented as a
% polytope Cx <= d_ub, where it also holds that Cx >= d_lb)

% read some values
n = dim(cZ);
nrEqcZ = length(cZ.b);
nrIneqP = size(C,2);

% same center
c_ = cZ.c;
% add zero-length generators
G = cZ.G;
G_ = [G, zeros(n,nrIneqP)];

% new constraint vector
A_ = cZ.A;
for i=1:nrIneqP
    con = C(:,i)';
    A_ = [A_, zeros(nrEqcZ+i-1,1);
          con*[G, zeros(n,i-1)], 0.5*(d_lb(i) - d_ub(i))];
end
% (note: c_ remains the same throughout the loop...)
b_ = [cZ.b; 0.5*(d_ub + d_lb) - C'*c_];

cZ_and = conZonotope(c_,G_,A_,b_);

end

% ------------------------------ END OF CODE ------------------------------
