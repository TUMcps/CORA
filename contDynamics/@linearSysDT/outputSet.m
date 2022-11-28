function Y = outputSet(obj,options,R)
% outputSet - calculates output set based on output equation given by
% 	 y = Cx + Du + k + v and sets for x (R) and u (options.U + uTrans)
%    and v (options.V)
%
% Syntax:  
%    Y = outputSet(obj,options,R)
%
% Inputs:
%    obj - linearSysDT object
%    options - options for the computation of reachable sets
%    R - reachable set of current step
%
% Outputs:
%    Y - time-point output set
%
% Example: 
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff, Mark Wetzlinger
% Written:      20-March-2020
% Last update:  16-November-2021 (MW, add sensor noise V)
% Last revision:---

%------------- BEGIN CODE --------------

% no output equation or output equation is y = x
if isempty(obj.C) || ...
        ( isscalar(obj.C) && obj.C == 1 && ~any(any(obj.D)) && ~any(obj.k) && ...
        ~any(center(options.V)) && isempty(generators(options.V)) )
    Y = R;
    return;
end

isD = false;
if any(any(obj.D))
    isD = true;
    U = options.U + options.uTrans;
end

if ~isfield(options,'saveOrder')
    % no additional reduction
    if isD
        Y = obj.C*R + obj.D * U + obj.k + options.V;
    else
        Y = obj.C*R + obj.k + options.V;
    end

else
    % reduction by saveOrder
    if isD
        Y = reduce(zonotope(obj.C*R) + obj.D * U + obj.k + options.V,...
            options.reductionTechnique,options.saveOrder);
    else
        Y = reduce(zonotope(obj.C*R) + obj.k + options.V,...
            options.reductionTechnique,options.saveOrder);
    end

end

%------------- END OF CODE --------------
