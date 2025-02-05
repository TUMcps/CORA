function timePoint = priv_reach_backward_EA_timepoint(linsys,params,options)
% priv_reach_backward_EA_timepoint - computes an inner approximation or an outer
%    approximation of the time-point maximal backward reachable set defined
%    as
%       R(-t) := { x_0 in R^n | forall u() exists w() :
%                               x(t; x_0, u(), w()) in Xend }
%    via [1, Sec. VI-A-1)]
%
% Syntax:
%    timePoint = priv_reach_backward_EA_timepoint(linsys,params,options)
%
% Inputs:
%    linsys - linearSys object
%    params - model parameters
%    options - reachability settings
%
% Outputs:
%    timePoint - struct of time-point solution
%
% References:
%    [1] M. Wetzlinger and M. Althoff, "Backward Reachability Analysis with
%        State Constraints for Linear Systems using Set Propagation", 2023.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: reachBackward

% Authors:       Mark Wetzlinger
% Written:       12-July-2023
% Last update:   28-November-2024 (MW, remove exponentialMatrix class)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% number of steps
steps = round(params.tFinal/options.timeStep);
% time step size (shorter variable name)
dt = options.timeStep;

% init helper class
linsys.taylor = taylorLinSys(linsys.A);

switch options.linAlg
    case 'inner:EA:timepoint'
        timePoint = aux_innerEAtimepoint(linsys,params.R0,params.U,...
            params.W,params.tFinal,dt,steps);
    case 'outer:EA:timepoint'
        timePoint = aux_outerEAtimepoint(linsys,params.R0,params.U,...
            params.W,params.tFinal,dt,steps);
end

end


% Auxiliary functions -----------------------------------------------------

function timePoint = aux_innerEAtimepoint(linsys,R0,U,W,tFinal,dt,steps)

% state dimension
n = dim(R0);

% compute outer- and inner-approximation of particular solutions [0,dt]
ZU_inner_0_dt = particularSolution_constant(linsys,U,dt,Inf);
ZW_outer_0_dt = particularSolution_timeVarying(linsys,W,dt,Inf);

% center and generator matrix of particular solution
if isa(ZW_outer_0_dt,'contSet')
    c_ZW = center(ZW_outer_0_dt);
    G_ZW = generators(ZW_outer_0_dt);
elseif isnumeric(ZW_outer_0_dt)
    c_ZW = ZW_outer_0_dt;
    G_ZW = zeros(linsys.nrOfDims,1);
end

% initialize support function values for particular solutions
l_0 = R0.A';

% initialize particular solutions
ZU_inner_0_t = zeros(n,1);
ZW_outer_0_t_sF = 0;

% pre-compute propagation matrices
eAtk = computeField(linsys.taylor,'eAdt',struct('timeStep',0));
eAdt = expm(linsys.A*dt);

for k=1:steps
    % back-propagate direction
    l = readFieldForTimeStep(linsys.taylor,'eAdt',dt*(k-1))' * l_0;

    % propagate particular solutions
    ZU_inner_0_t = ZU_inner_0_t + ...
        readFieldForTimeStep(linsys.taylor,'eAdt',dt*(k-1)) * ZU_inner_0_dt;
    ZW_outer_0_t_sF = ZW_outer_0_t_sF + (l'*c_ZW + sum(abs(l' * G_ZW),2));

    % propagate propagation matrix
    eAtk = eAtk*eAdt;
    insertFieldTimeStep(linsys.taylor,'eAdt',eAtk,dt*k);
end

S_inner_0_t = polytope(R0.A,R0.b - ZW_outer_0_t_sF);

if isemptyobject(S_inner_0_t)
    timePoint.set{1} = interval.empty(n);
else
    timePoint.set{1} = expm(-linsys.A*tFinal) ...
            * (conZonotope(S_inner_0_t) + -1*ZU_inner_0_t);
end
timePoint.time{1} = tFinal;

end

function timePoint = aux_outerEAtimepoint(linsys,R0,U,W,tFinal,dt,steps)

% read out dimension
n = dim(R0);

% compute outer- and inner-approximation of particular solutions [0,dt]
ZU_outer_0_dt = particularSolution_timeVarying(linsys,U,dt,Inf);
ZW_inner_0_dt = particularSolution_constant(linsys,W,dt,Inf);

% center and generator matrix of particular solution
if isa(ZW_inner_0_dt,'contSet')
    c_ZW = center(ZW_inner_0_dt);
    G_ZW = generators(ZW_inner_0_dt);
elseif isnumeric(ZW_inner_0_dt)
    c_ZW = ZW_inner_0_dt;
    G_ZW = zeros(linsys.nrOfDims,1);
end

% initialize support function values for particular solutions
l_0 = R0.A';
ZW_inner_0_t_sF = 0;
ZU_outer_0_t = zeros(n,1);

% pre-compute propagation matrices
eAtk = computeField(linsys.taylor,'eAdt',struct('timeStep',0));
eAdt = expm(linsys.A*dt);

% compute Minkowski difference Xend - PU/PW for all steps
for k=1:steps
    % propagate directions
    l = readFieldForTimeStep(linsys.taylor,'eAdt',dt*(k-1))' * l_0;

    % propagate particular solutions
    ZU_outer_0_t = ZU_outer_0_t + ...
        readFieldForTimeStep(linsys.taylor,'eAdt',dt*(k-1)) * ZU_outer_0_dt;
    ZW_inner_0_t_sF = ZW_inner_0_t_sF + (l'*c_ZW + sum(abs(l' * G_ZW),2));

    % propagate propagation matrix
    eAtk = eAtk*eAdt;
    insertFieldTimeStep(linsys.taylor,'eAdt',eAtk,dt*k);
end

S_outer_0_t = polytope(R0.A,R0.b - ZW_inner_0_t_sF);

% Minkowski addition with (-1) times other uncertainty set; then linear map
if isemptyobject(S_outer_0_t)
    timePoint.set{1} = interval.empty(n);
else
    timePoint.set{1} = expm(-linsys.A*tFinal) ...
        * (conZonotope(S_outer_0_t) + -1*ZU_outer_0_t);
end
timePoint.time{1} = tFinal;

end

% ------------------------------ END OF CODE ------------------------------
