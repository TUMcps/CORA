function [H,Zdelta,errorStat,T,ind3,Zdelta3] = precompStatError_adaptive(obj,options,Rdelta)
% precompStatError_adaptive - precompute the second order static error
%    along with Hessian matrix, with adaptive zonotope order reduction
%
% Syntax:
%    [H,Zdelta,errorStat,T,ind3,Zdelta3] = precompStatError_adaptive(obj,options,Rdelta)
%
% Inputs:
%    obj - nonlinear system object
%    Rdelta - shifted reachable set at the beginning of the time step
%    options - options struct
%
% Outputs:
%    H - hessian matrix
%    Zdelta - zonotope over-approximating the reachable set at the
%             beginning of the time step extended by the input set 
%    errorStat - static linearization error
%    T - third-order tensor
%    ind3 - indices at which the third-order tensor is not zero
%    Zdelta3 - set Zdelta reduced to the zonotope order for the evaluation
%              of the third-order tensor
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% References: 
%   [1] M. Althoff et al. "Reachability Analysis of Nonlinear Systems with 
%       Uncertain Parameters using Conservative Linearization"
%   [2] M. Althoff et al. "Reachability analysis of nonlinear systems using 
%       conservative polynomialization and non-convex sets"
%
% See also: linReach, linError_*, preCompStatError

% Authors:       Mark Wetzlinger
% Written:       ---
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% precompute the second order static error along with hessian matrix
obj = setHessian(obj,'standard');

% initialize output arguments
T = []; ind3 = []; Zdelta3 = [];

% reduce the reachable set for the initial time point
Rred = reduce(Rdelta,'adaptive',sqrt(options.redFactor));

% over-approximation of the reachable set at the initial time point
Rdelta = reduce(zonotope(Rdelta),'adaptive',sqrt(options.redFactor));

% extend the sets by the input sets
Z = cartProd(Rred,options.U);
Zdelta = cartProd(Rdelta,options.U);

% calculate the hessian tensor
H = obj.hessian(obj.linError.p.x, obj.linError.p.u);

% calculate the quadratic map == static second-order error
errorSecOrdStat = 0.5*quadMap(Z,H);

% reduce the order of the set of static errors
errorStat = reduce(errorSecOrdStat,'adaptive',sqrt(options.redFactor));

% ------------------------------ END OF CODE ------------------------------
