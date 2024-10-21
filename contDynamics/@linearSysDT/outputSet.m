function Y = outputSet(linsysDT,R,params,options)
% outputSet - calculates output set based on output equation given by
% 	 y = Cx + Du + k + v and sets for x (R) and u (params.U + params.uTrans)
%    and v (options.V)
%
% Syntax:
%    Y = outputSet(linsysDT,R,params,options)
%
% Inputs:
%    linsysDT - linearSysDT object
%    R - reachable set of current step
%    params - model parameters
%    options - options for the computation of reachable sets
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

% Authors:       Matthias Althoff, Mark Wetzlinger
% Written:       20-March-2020
% Last update:   16-November-2021 (MW, add sensor noise V)
%                07-December-2022 (MW, allow to skip output set)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% skip computation of output set
if ~options.compOutputSet
    Y = R; return
end

% no output equation or output equation is y = x
if isempty(linsysDT.C) || ...
        ( isscalar(linsysDT.C) && linsysDT.C == 1 && ~any(any(linsysDT.D)) ...
        && ~any(linsysDT.k) && representsa_(params.V,'origin',eps) )
    Y = R;
    return;
end

isD = false;
if iscell(linsysDT.D) || any(any(linsysDT.D))
    isD = true;
    U = params.U + params.uTrans;
end

if ~isfield(options,'saveOrder')
    % no additional reduction
    if iscell(linsysDT.C)
        Y = linsysDT.C{options.i}*R + linsysDT.D{options.i} * U + linsysDT.k + params.V;
    else
        if isD
            Y = linsysDT.C*R + linsysDT.D * U + linsysDT.k + params.V;
        else
            Y = linsysDT.C*R + linsysDT.k + params.V;
        end
    end

else
    % reduction by saveOrder
    if iscell(linsysDT.C)
        Y = reduce(zonotope(linsysDT.C{options.i}*R) + linsysDT.D{options.i} * U + linsysDT.k + params.V,...
                options.reductionTechnique,options.saveOrder);
    else
        if isD
            Y = reduce(zonotope(linsysDT.C*R) + linsysDT.D * U + linsysDT.k + params.V,...
                options.reductionTechnique,options.saveOrder);
        else
            Y = reduce(zonotope(linsysDT.C*R) + linsysDT.k + params.V,...
                options.reductionTechnique,options.saveOrder);
        end
    end

end

% ------------------------------ END OF CODE ------------------------------
