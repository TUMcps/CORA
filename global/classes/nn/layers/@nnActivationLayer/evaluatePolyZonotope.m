function [c, G, Grest, expMat, id, id_, ind, ind_] = evaluatePolyZonotope(obj, c, G, Grest, expMat, id, id_, ind, ind_, evParams)
% evaluatePolyZonotope - evaluates the activation layer on a polyZonotope
%
% Syntax:
%    [c, G, Grest, expMat, id, id_, ind, ind_] = evaluatePolyZonotope(obj, c, G, Grest, expMat, id, id_, ind, ind_, evParams)
%
% Inputs:
%    c, G, Grest, expMat, id, id_, ind, ind_ - parameters of polyZonotope
%    evParams - parameter for NN evaluation
%
% Outputs:
%    updated [c, G, Grest, expMat, id, id_, ind, ind_]
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: nnLayer

% Author:        Niklas Kochdumper, Tobias Ladner
% Written:       17-September-2021
% Last update:   ---
% Last revision: 05-April-2022 (TL)

%------------- BEGIN CODE --------------

% LINEAR APPROXIMATION
if strcmp(evParams.polynomial_approx, "lin")
    [c, G, Grest, expMat, id, id_, ind, ind_] = evaluatePolyZonotopeLin(obj, c, G, Grest, expMat, id, id_, ind, ind_, evParams);

    % QUADRADIC APPROXIMATION
elseif strcmp(evParams.polynomial_approx, "quad")
    [c, G, Grest, expMat, id, id_, ind, ind_] = evaluatePolyZonotopeQuad(obj, c, G, Grest, expMat, id, id_, ind, ind_, evParams);

    % CUBIC APPROXIMATION
elseif strcmp(evParams.polynomial_approx, "cub")
    [c, G, Grest, expMat, id, id_, ind, ind_] = evaluatePolyZonotopeCub(obj, c, G, Grest, expMat, id, id_, ind, ind_, evParams);
end

%------------- END OF CODE --------------