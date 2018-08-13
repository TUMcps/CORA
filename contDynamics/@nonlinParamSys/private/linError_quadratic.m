function [error,errorInt] = linError_quadratic(obj,options,R)
% linError - computes the linearization error
%
% Syntax:  
%    [obj] = linError(obj,options)
%
% Inputs:
%    obj - nonlinear system object
%    options - options struct
%    R - actual reachable set
%
% Outputs:
%    obj - nonlinear system object
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: 

% Author:       Matthias Althoff
% Written:      29-October-2007 
% Last update:  22-January-2008
%               02-February-2010
%               25-July-2016 (intervalhull replaced by interval)
%               07-July-2017
% Last revision: ---

%------------- BEGIN CODE --------------

%compute interval of reachable set
totalInt=interval(R) + obj.linError.p.x;

%compute intervals of input
inputInt=interval(options.U) + options.uTrans;

%obtain hessian tensor
H = obj.hessian(totalInt,inputInt,options.paramInt);

%compute zonotope of state and input
Rred = reduce(R,options.reductionTechnique,options.errorOrder);
Z=cartesianProduct(Rred,options.U);

%compute input due to second order evaluation
error=0.5*quadraticMultiplication_interval(Z,H);


error = reduce(error,options.reductionTechnique,options.zonotopeOrder);

errorIHabs = abs(interval(error));
errorInt = supremum(errorIHabs);


%------------- END OF CODE --------------