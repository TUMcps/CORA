function timePoint = reach_backward_AE_timepoint(linsys,params,options)
% reach_backward_AE_timepoint - computes an inner approximation or an
%    outer approximation of the time-point minimal backward reachable set
%    defined as
%       R(-t) := { x_0 in R^n | exists u() forall w() :
%                               x(t; x_0, u(), w()) in Xend }
%    and computed via [1, Sec. V-A]
%
% Syntax:
%    timePoint = reach_backward_AE_timepoint(linsys,params,options)
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

% Authors:        Mark Wetzlinger
% Written:        12-July-2023
% Last update:    ---
% Last revision:  ---

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
    case 'outer:AE:timepoint'
        timePoint = aux_outerAEtimepoint(expmat,params.R0,params.U,...
            params.W,params.tFinal,dt,steps);
    case 'inner:AE:timepoint'
        timePoint = aux_innerAEtimepoint(expmat,params.R0,params.U,...
            params.W,params.tFinal,dt,steps);
end

end


% Auxiliary functions -----------------------------------------------------

function timePoint = aux_outerAEtimepoint(expmat,R0,U,W,tFinal,dt,steps)

% compute inner- and outer-approximations of particular solution [0,dt]
% (save powers of A in the process...)
[ZW_outer_0_dt,expmat] = oneStepPartSol(W,expmat,dt,'outer');
[ZU_inner_0_dt,expmat] = oneStepPartSol(U,expmat,dt,'inner');

% centers and generator matrices of particular solutions
c_ZW = -center(ZW_outer_0_dt);
G_ZW = -generators(ZW_outer_0_dt);
c_ZU = center(ZU_inner_0_dt);
G_ZU = generators(ZU_inner_0_dt);

% initialize support function values for particular solutions
ZW_outer_0_t_sF = 0;
ZU_inner_0_t_sF = 0;
l_0 = R0.A';

% propagation of direction and support function values
for k=1:steps
    % back-propagate direction
    expmat = next_eAtk(expmat,k,'transpose');
    l = expmat.eAtk_T{k} * l_0;

    % propagate support function values of particular solutions
    ZW_outer_0_t_sF = ZW_outer_0_t_sF + ...
        l'*c_ZW + sum(abs(l' * G_ZW),2);
    ZU_inner_0_t_sF = ZU_inner_0_t_sF + ...
        l'*c_ZU + sum(abs(l' * G_ZU),2);
end

% compute time-point solution
timePoint.set{1} = expm(-expmat.A*tFinal) * ...
    polytope(R0.A,R0.b - ZU_inner_0_t_sF + ZW_outer_0_t_sF);
timePoint.time{1} = tFinal;

end

function timePoint = aux_innerAEtimepoint(expmat,R0,U,W,tFinal,dt,steps)

% compute inner- and outer-approximations of particular solution [0,dt]
[ZW_inner_0_dt,expmat] = oneStepPartSol(W,expmat,dt,'inner');
[ZU_outer_0_dt,expmat] = oneStepPartSol(U,expmat,dt,'outer');

% center and generator matrix of particular solution
c_ZU = center(ZU_outer_0_dt);
G_ZU = generators(ZU_outer_0_dt);

% initialize support function values for particular solutions
ZU_outer_0_t_sF = 0;
l_0 = R0.A';

% init inner-approximation of particular solution
ZW_inner_0_t = zonotope(zeros(dim(R0),1));

for k=1:steps
    % back-propagate direction
    expmat = next_eAtk(expmat,k,'transpose');
    l = expmat.eAtk_T{k} * l_0;

    % propagate support function value
    ZU_outer_0_t_sF = ZU_outer_0_t_sF + (l'*c_ZU + sum(abs(l' * G_ZU),2));
    
    % propagate particular solution
    expmat = next_eAtk(expmat,k);
    ZW_inner_0_t = ZW_inner_0_t + expmat.eAtk{k} * ZW_inner_0_dt;        
end

% compute Minkowski difference Xend - PU
temp_inner_0_t = polytope(R0.A,R0.b - ZU_outer_0_t_sF);

% Minkowski addition with (-1) times other uncertainty set; then linear map
timePoint.set{1} = expm(-expmat.A*tFinal) ...
    * (conZonotope(temp_inner_0_t) + -1*ZW_inner_0_t);
timePoint.time{1} = tFinal;

end

% ------------------------------ END OF CODE ------------------------------
