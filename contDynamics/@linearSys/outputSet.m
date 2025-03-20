function Y = outputSet(linsys,R,params,options)
% outputSet - calculates output set based on output equation given by
%    y = Cx + Du + k + Fv and sets for x (R) and u (options.U + options.uTrans)
%
% Syntax:
%    Y = outputSet(linsys,R,params,options)
%
% Inputs:
%    linsys - linearSys object
%    R - reachable set (either time point [i] or time interval [i,i+1])
%    params - model parameters
%    options - options for the computation of reachable sets
%
% Outputs:
%    Y - output set (either time point [i] or time interval [i,i+1])
%
% Example:
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
if isempty(linsys.C) || ...
        ( isscalar(linsys.C) && linsys.C == 1 && ~any(any(linsys.D)) ...
        && ~any(linsys.k) && (~any(any(linsys.F)) || representsa_(params.V,'origin',eps)) )
    Y = R;
    return;
end


% do we consider inputs in the output?
isD = false;
if any(any(linsys.D))
    isD = true;
    U = params.U + params.uTrans;
end


if ~isfield(options,'saveOrder')
    % output equations without order reduction
    if isD
        Y = linsys.C*R + linsys.D * U + linsys.k + linsys.F * params.V;
    else
        Y = linsys.C*R + linsys.k + linsys.F * params.V;
    end

else

    % reduction by saveOrder
    if isD
        Y = reduce(zonotope(linsys.C*R) + linsys.D * U + linsys.k + linsys.F * params.V,...
            options.reductionTechnique,options.saveOrder);
    else
        Y = reduce(zonotope(linsys.C*R) + linsys.k + linsys.F * params.V,...
            options.reductionTechnique,options.saveOrder);
    end

end

% ------------------------------ END OF CODE ------------------------------
