function [error,errorInt] = linError_mixed_noInt(obj,options,R)
% linError_mixed_noInt - computes the linearization error without use of
% interval arithmatic
%
% Syntax:  
%    [error,errorInt] = linError_mixed_noInt(obj,options,R)
%
% Inputs:
%    obj - nonlinear system object
%    options - options struct
%    R - actual reachable set
%
% Outputs:
%    error - error represented by a zonotope
%    errorInt - multidimensional interval of error
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: 

% Author:       Matthias Althoff
% Written:      11-July-2012
% Last update:  25-July-2016 (intervalhull replaced by interval)
%               07-July-2017
% Last revision:---

%------------- BEGIN CODE --------------


%compute interval of reachable set
IH_x = interval(R);
totalInt_x = interval(IH_x) + obj.linError.p.x;

%compute intervals of input
IH_u = interval(options.U);
totalInt_u = interval(IH_u) + obj.linError.p.u;

%compute zonotope of state and input
Rred = reduce(R,options.reductionTechnique,options.errorOrder);
Z=cartesianProduct(Rred,options.U);

%obtain hessian tensor
if isfield(options,'lagrangeRem') && isfield(options.lagrangeRem,'method') && ...
   ~strcmp(options.lagrangeRem.method,'interval')

    % create taylor models or zoo-objects
    [objX,objU] = initRangeBoundingObjects(totalInt_x,totalInt_u,options);

    % evaluate the hessian tensor 
    H = obj.hessian(objX,objU,options.paramInt);
else
    H = obj.hessian(totalInt_x, totalInt_u, options.paramInt);
end

%obtain intervals and combined interval z
dx = interval(IH_x);
du = interval(IH_u);
dz = [dx; du];

%obtain absolute values
dz_abs = max(abs(infimum(dz)), abs(supremum(dz)));

%separate evaluation
for i=1:length(H)
    H_mid{i} = sparse(mid(H{i}));
    H_rad{i} = sparse(rad(H{i}));
end

error_mid = 0.5*quadraticMultiplication(Z, H_mid);

%interval evaluation
for i=1:length(H)
    error_rad(i,1) = 0.5*dz_abs'*H_rad{i}*dz_abs;
end

%combine results
error_rad_zono = zonotope(interval(-error_rad, error_rad));
error = error_mid + error_rad_zono;

if isa(error, 'zonotope')
    error = reduce(error,options.reductionTechnique,options.zonotopeOrder);
else
    error = reduce(error,options.reductionTechnique,options.intermediateOrder);
end

errorIHabs = abs(interval(error));
errorInt = supremum(errorIHabs);

% H2 = hessianTensor(mid(totalInt),mid(inputInt));
% error2=0.5*quadraticMultiplication(Z,H2);



%------------- END OF CODE --------------