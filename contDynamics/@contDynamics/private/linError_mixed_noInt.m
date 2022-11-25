function [errorZon,errorInt] = linError_mixed_noInt(obj,options,R)
% linError_mixed_noInt - computes the linearization error without use of
% interval arithmetic (see Theorem 1 in [1])
%
% Syntax:  
%    [error,errorInt] = linError_mixed_noInt(obj,options,R)
%
% Inputs:
%    obj - nonlinear system object
%    options - options struct
%    R - reachable set including admissible error
%
% Outputs:
%    errorZon - error represented by a zonotope
%    errorInt - multidimensional interval of error
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
%   [1] M. Althoff and B. Krogh
%       "Reachability Analysis of Nonlinear Differential-Algebraic Systems"

% Author:       Matthias Althoff
% Written:      11-July-2012
% Last update:  25-July-2016 (intervalhull replaced by interval)
%               21-April-2020 (added reference, simplification)
% Last revision:---

%------------- BEGIN CODE --------------


%compute interval of reachable set
IH_x = interval(R);
totalInt_x = IH_x + obj.linError.p.x;

%compute intervals of input
IH_u = interval(options.U);
totalInt_u = IH_u + obj.linError.p.u;

%compute zonotope of state and input
Rred = reduce(zonotope(R),options.reductionTechnique,options.errorOrder);
Z = cartProd(Rred,options.U);

%obtain hessian tensor
if isfield(options,'lagrangeRem') && isfield(options.lagrangeRem,'method') && ...
   ~strcmp(options.lagrangeRem.method,'interval')

    % create taylor models or zoo-objects
    [objX,objU] = initRangeBoundingObjects(totalInt_x,totalInt_u,options);

    % evaluate the hessian tensor
    if isa(obj,'nonlinParamSys')
        H = obj.hessian(objX,objU,options.paramInt);
    else
        H = obj.hessian(objX,objU);
    end
else
    if isa(obj,'nonlinParamSys')
        H = obj.hessian(totalInt_x,totalInt_u,options.paramInt);
    else
        H = obj.hessian(totalInt_x,totalInt_u);
    end
end

%obtain combined interval z and absolute values
dz = [IH_x; IH_u];
dz_abs = max(abs(infimum(dz)),abs(supremum(dz)));

%separate evaluation
H_mid = cell(obj.dim,1);
H_rad = cell(obj.dim,1);
for i=1:obj.dim
    H_mid{i} = sparse(center(H{i}));
    H_rad{i} = sparse(rad(H{i}));
end

error_mid = 0.5*quadMap(Z,H_mid);

%interval evaluation
error_rad = zeros(obj.dim,1);
for i=1:obj.dim
    error_rad(i,1) = 0.5*dz_abs'*H_rad{i}*dz_abs;
end

%combine results
error_rad_zono = zonotope(interval(-error_rad,error_rad));
errorZon = error_mid + error_rad_zono;
errorZon = reduce(errorZon,options.reductionTechnique,options.intermediateOrder);

errorInt = supremum(abs(interval(errorZon)));

%------------- END OF CODE --------------