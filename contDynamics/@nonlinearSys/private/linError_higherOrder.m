function [Rerror,Error] = linError_higherOrder(obj,R,options)
% linError_higherOrder - computes the linearization error by using higher
%                        order tensors
%
% Syntax:  
%    [Rerror,Error] = linError_higherOrder(obj,R,options)
%
% Inputs:
%    obj - nonlinear system object
%    R - reachable set of the current time interval
%    options - options struct
%
% Outputs:
%    Rerror - set of linearization errors (class: zonotope)
%    Error - upper boundary for the absolute error
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: 

% Author:       Niklas Kochdumper
% Written:      02-March-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% compute interval of reachable set
dx = interval(R);
totalInt_x = dx + obj.linError.p.x;

% compute intervals of input
du = interval(options.U);
totalInt_u = du + obj.linError.p.u;

% obtain intervals and combined interval z
dz = [dx; du];

% reduce zonotope
R = reduce(R,options.reductionTechnique,options.errorOrder);

% combined zonope (states + input)
Z = cartesianProduct(R,options.U);

% calculate hessian matrix
H = obj.hessian(obj.linError.p.x,obj.linError.p.u);

% calculate third-order tensor
if isfield(options,'lagrangeRem') && isfield(options.lagrangeRem,'method') && ...
    ~strcmp(options.lagrangeRem.method,'interval')

    % create taylor models or zoo-objects
    [objX,objU] = initRangeBoundingObjects(totalInt_x,totalInt_u,options);

    % evaluate third order tensor 
    T = obj.thirdOrderTensor(objX, objU);

else
    T = obj.thirdOrderTensor(totalInt_x, totalInt_u);
end

% error second-order
errorSec = 0.5 * quadraticMultiplication(Z,H);

% calculate the Lagrange remainder term
for i=1:length(T(:,1))
    error_sum = interval(0,0);
    for j=1:length(T(1,:))
        error_tmp(i,j) = dz.'*T{i,j}*dz;
        error_sum = error_sum + error_tmp(i,j) * dz(j);
    end
    errorLagr(i,1) = 1/6*error_sum;
end

% overall linearization error
Rerror = errorSec + zonotope(errorLagr);

Error = abs(interval(Rerror));
Error = supremum(Error);

%------------- END OF CODE --------------