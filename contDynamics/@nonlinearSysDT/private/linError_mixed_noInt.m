function error = linError_mixed_noInt(obj, options, R)
% linError_mixed_noInt - computes the linearization error
%
% Syntax:  
%    error = linError_mixed_noInt(obj,options,R)
%
% Inputs:
%    obj - nonlinearSysDT system object
%    options - options struct
%    R - initial set
%
% Outputs:
%    error - zonotope overapproximating the linearization error
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: 

% Author:       Matthias Althoff, Niklas Kochdumper
% Written:      21-August-2012
% Last update:  25-July-2016 (intervalhull replaced by interval)
%               29-January-2018 (NK)
% Last revision:---

%------------- BEGIN CODE --------------

%obtain intervals and combined interval z
dx = interval(R);
du = interval(options.U);
dz = [dx; du];

%compute interval of reachable set
totalInt_x = dx + obj.linError.p.x;

%compute intervals of input
totalInt_u = du + obj.linError.p.u;

%compute zonotope of state and input
Rred = reduce(R,'girard',options.errorOrder);
Z=cartesianProduct(Rred,options.U);

%obtain hessian tensor
if isfield(options,'lagrangeRem') && isfield(options.lagrangeRem,'method') && ...
   ~strcmp(options.lagrangeRem.method,'interval')

    % create taylor models or zoo-objects
    [objX,objU] = initRangeBoundingObjects(totalInt_x,totalInt_u,options);

    % evaluate the hessian tensor 
    H = obj.hessian(objX,objU,options.timeStep);
else
    H = obj.hessian(totalInt_x, totalInt_u, options.timeStep);
end

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

error = reduce(error,'girard',options.zonotopeOrder);

%------------- END OF CODE --------------