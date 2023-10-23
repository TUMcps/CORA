function Y = outputSet(obj,options,R)
% outputSet - calculates output set based on output equation given by
%    y = Cx + Du + k + v and sets for x (R) and u (options.U + options.uTrans)
%
% Syntax:
%    Y = outputSet(obj,options,R)
%
% Inputs:
%    obj - linearSys object
%    options - options for the computation of reachable sets
%    R - reachable set (either time point [i] or time interval [i,i+1])
%
% Outputs:
%    Y - output set (either time point [i] or time interval [i,i+1])
%
% Example:
%    -
%
% References:
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       12-August-2019
% Last update:   20-August-2019
%                16-November-2021 (MW, add sensor noise V)
%                19-November-2021 (MW, shift index of time-point solution)
%                19-November-2022 (MW, remove double computation)
%                07-December-2022 (MW, allow to skip output set)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% skip computation of output set
if ~options.compOutputSet
    Y = R; return
end

% output equation is not provided or y = x
if isempty(obj.C) || ...
        ( isscalar(obj.C) && obj.C == 1 && ~any(any(obj.D)) ...
        && ~any(obj.k) && representsa_(options.V,'origin',eps) )
    Y = R;
    return;
end

isD = false;
if any(any(obj.D))
    isD = true;
    U = options.U + options.uTrans;
end


if ~isfield(options,'saveOrder')
    
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

% ------------------------------ END OF CODE ------------------------------
