function timePoint = reach_backward_EA_timepoint(linsys,params,options)
% reach_backward_EA_timepoint - computes an inner approximation or an outer
%    approximation of the time-point maximal backward reachable set defined
%    as
%       R(-t) := { x_0 in R^n | forall u() exists w() :
%                               x(t; x_0, u(), w()) in Xend }
%    via [1, Sec. VI-A-1)]
%
% Syntax:
%    timePoint = reach_backward_EA_timepoint(linsys,params,options)
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
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% number of steps
steps = round(params.tFinal/options.timeStep);
% time step size (shorter variable name)
dt = options.timeStep;

% initialize exponential matrix
expmat = exponentialMatrix(linsys.A);
% time step size known, compute propagation matrix over one time step
expmat = init_eADeltat(expmat,dt);
% init cell length for propagation matrices
expmat = init_eAtk(expmat,steps);

switch options.linAlg
    case 'inner:EA:timepoint'
        timePoint = aux_innerEAtimepoint(expmat,params.R0,params.U,...
            params.W,params.tFinal,dt,steps);
    case 'outer:EA:timepoint'
        timePoint = aux_outerEAtimepoint(expmat,params.R0,params.U,...
            params.W,params.tFinal,dt,steps);
end

end


% Auxiliary functions -----------------------------------------------------

function timePoint = aux_innerEAtimepoint(expmat,R0,U,W,tFinal,dt,steps)

% state dimension
n = dim(R0);

% compute outer- and inner-approximation of particular solutions [0,dt]
[ZU_inner_0_dt,expmat] = oneStepPartSol(U,expmat,dt,'inner');
[ZW_outer_0_dt,expmat] = oneStepPartSol(W,expmat,dt,'outer');

% center and generator matrix of particular solution
c_ZW = center(ZW_outer_0_dt);
G_ZW = generators(ZW_outer_0_dt);

% initialize support function values for particular solutions
l_0 = R0.A';

% initialize particular solutiosn
ZU_inner_0_t = zeros(n,1);
ZW_outer_0_t_sF = 0;

for k=1:steps
    % propagate directions
    expmat = next_eAtk(expmat,k,'transpose');
    l = expmat.eAtk_T{k} * l_0;

    % propagate particular solutions
    expmat = next_eAtk(expmat,k);
    ZU_inner_0_t = ZU_inner_0_t + expmat.eAtk{k} * ZU_inner_0_dt;
    ZW_outer_0_t_sF = ZW_outer_0_t_sF + (l'*c_ZW + sum(abs(l' * G_ZW),2));
end

temp_inner_0_t = polytope(R0.A,R0.b - ZW_outer_0_t_sF);

if isemptyobject(temp_inner_0_t)
    timePoint.set{1} = interval(zeros(n,0),zeros(n,0));
else
    timePoint.set{1} = expm(-expmat.A*tFinal) ...
            * (conZonotope(temp_inner_0_t) + -1*ZU_inner_0_t);
end
timePoint.time{1} = tFinal;

end

function timePoint = aux_outerEAtimepoint(expmat,R0,U,W,tFinal,dt,steps)

% read out dimension
n = dim(R0);

% compute outer- and inner-approximation of particular solutions [0,dt]
[ZU_outer_0_dt,expmat] = oneStepPartSol(U,expmat,dt,'outer');
[ZW_inner_0_dt,expmat] = oneStepPartSol(W,expmat,dt,'inner');

% center and generator matrix of particular solution
c_ZW = center(ZW_inner_0_dt);
G_ZW = generators(ZW_inner_0_dt);

% initialize support function values for particular solutions
l_0 = R0.A';
ZW_inner_0_t_sF = 0;
ZU_outer_0_t = zeros(n,1);

% compute Minkowski difference Xend - PU/PW for all steps
for k=1:steps
    % propagate directions
    expmat = next_eAtk(expmat,k,'transpose');
    l = expmat.eAtk_T{k} * l_0;

    % propagate particular solutions
    expmat = next_eAtk(expmat,k);
    ZU_outer_0_t = ZU_outer_0_t + expmat.eAtk{k} * ZU_outer_0_dt;
    ZW_inner_0_t_sF = ZW_inner_0_t_sF + (l'*c_ZW + sum(abs(l' * G_ZW),2));
end

temp_outer_0_t = polytope(R0.A,R0.b - ZW_inner_0_t_sF);

% Minkowski addition with (-1) times other uncertainty set; then linear map
if isemptyobject(temp_outer_0_t)
    timePoint.set{1} = interval.empty(n);
else
    timePoint.set{1} = expm(-expmat.A*tFinal) ...
        * (conZonotope(temp_outer_0_t) + -1*ZU_outer_0_t);
end
timePoint.time{1} = tFinal;

end

% ------------------------------ END OF CODE ------------------------------
