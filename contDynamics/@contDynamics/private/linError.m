function [error] = linError(obj,options,R)
% linError - computes the linearization error
%
% Syntax:  
%    [error] = linError(obj,options)
%
% Inputs:
%    obj - nonlinearSys or nonlinParamSys object
%    options - options struct
%    R - actual reachable set
%
% Outputs:
%    error - linearization error
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
%
% See also: linReach

% Author:       Matthias Althoff, Niklas Kochdumper
% Written:      29-October-2007 
% Last update:  22-January-2008
%               02-February-2010
%               25-July-2016 (intervalhull replaced by interval)
%               12-November-2018 (NK: changed method for remainder
%                                 over-approximation)
% Last revision: ---

%------------- BEGIN CODE --------------

% compute interval of reachable set
IHx = interval(R);
% compute intervals of total reachable set
totalInt = IHx + obj.linError.p.x;

% compute intervals of input
inputInt = interval(options.U) + options.uTrans;
% translate intervals by linearization point
IHu = inputInt + (-obj.linError.p.u);

% obtain maximum absolute values within IH, IHinput
dx = max(abs(infimum(IHx)),abs(supremum(IHx)));
du = max(abs(infimum(IHu)),abs(supremum(IHu)));

% evaluate the hessian matrix with the selected range-bounding technique
if isfield(options,'lagrangeRem') && isfield(options.lagrangeRem,'method') && ...
   ~strcmp(options.lagrangeRem.method,'interval')

    % create taylor models or zoo-objects
    [objX,objU] = initRangeBoundingObjects(totalInt,inputInt,options);

    % evaluate the Lagrange remainder 
    if isa(obj,'nonlinParamSys')
        H = obj.hessian(objX,objU,options.paramInt);
    else
        H = obj.hessian(objX,objU);
    end
else
    if isa(obj,'nonlinParamSys')
        H = obj.hessian(totalInt,inputInt,options.paramInt);
    else
        H = obj.hessian(totalInt,inputInt);
    end
end

% over-approximate the Lagrange remainder acc. to Proposition 1 in [1]
error = zeros(length(H),1);
dz = [dx;du];

for i = 1:length(H)
    H_ = abs(H{i});
    H_ = max(infimum(H_),supremum(H_));
    error(i) = 0.5 * dz' * H_ * dz;
end

%------------- END OF CODE --------------