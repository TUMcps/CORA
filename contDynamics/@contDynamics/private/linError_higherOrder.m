function [Rerror,Error] = linError_higherOrder(obj,R,options)
% linError_higherOrder - computes the linearization error by using higher
%    order tensors (extension of the error computation in [1])
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
%
% References: 
%   [1] M. Althoff, O. Stursberg, and M. Buss
%       "Reachability Analysis of Nonlinear Systems with uncertain
%           Parameters using Conservative Linearization"

% Author:       Niklas Kochdumper
% Written:      02-March-2018
% Last update:  21-April-2020 (skip empty T{i,j})
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
Rred = reduce(R,options.reductionTechnique,options.errorOrder);

% combined zonotope (states + input)
Z = cartProd(Rred,options.U);

% calculate hessian matrix
if isa(obj,'nonlinParamSys')
    H = obj.hessian(obj.linError.p.x,obj.linError.p.u,options.paramInt);
else
    H = obj.hessian(obj.linError.p.x,obj.linError.p.u);
end

% evaluate third-order tensor 
if isfield(options,'lagrangeRem') && isfield(options.lagrangeRem,'method') && ...
    ~strcmp(options.lagrangeRem.method,'interval')

    % create taylor models or zoo objects
    [objX,objU] = initRangeBoundingObjects(totalInt_x,totalInt_u,options);

    if isa(obj,'nonlinParamSys')
        [T,ind] = obj.thirdOrderTensor(objX, objU, options.paramInt);
    else
        [T,ind] = obj.thirdOrderTensor(objX, objU);
    end

else
    if isa(obj,'nonlinParamSys')
        [T,ind] = obj.thirdOrderTensor(totalInt_x, totalInt_u, options.paramInt);
    else
        [T,ind] = obj.thirdOrderTensor(totalInt_x, totalInt_u);
    end
end

% second-order error
errorSec = 0.5 * quadMap(Z,H);

% calculate the Lagrange remainder term (third order)
if isempty(cell2mat(ind))
    % empty zonotope if all entries in T are empty
    errorLagr = zonotope(zeros(obj.dim,1));
else
    % skip tensors with all-zero entries using ind from tensor creation
    errorLagr = interval(zeros(obj.dim,1),zeros(obj.dim,1));
    for i=1:length(ind)
        error_sum = interval(0,0);
        for j=1:length(ind{i})
            error_tmp = dz.'*T{i,j}*dz;
            error_sum = error_sum + error_tmp * dz(j);
        end
        errorLagr(i,1) = 1/6*error_sum;
    end
    errorLagr = zonotope(errorLagr);
end

% overall linearization error
Rerror = errorSec + errorLagr;
Rerror = reduce(Rerror,options.reductionTechnique,options.intermediateOrder);

Error = abs(interval(Rerror));
Error = supremum(Error);

%------------- END OF CODE --------------