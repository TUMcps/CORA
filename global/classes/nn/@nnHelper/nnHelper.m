classdef nnHelper
% nnHelper - Helper class for neuralNetwork for easier usage of common
% functionalities
% Note: Many functions come from the neuralNetworkOld class.
% Some functions can probably be removed (unused legacy)
%
% Syntax:
%    only has static functions.
%
% Inputs:
%    -
%
% Outputs:
%    -
%
% Example:
%    nnHelper.myHelperfunction(...)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none

% Author:       Tobias Ladner
% Written:      28-March-2022
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

methods (Static)
    [l, u] = compBoundsPolyZono(c, G, Grest, E, ind, ind_, approx)
    [c, G, C, d, l, u] = conversionConZonoStarSet(cZ)
    cZ = conversionStarSetConZono(c, G, C, d, l_, u_)
    [c, G, Grest] = cubApproxPolyZono(c, G, Grest, c_a, c_b, c_c, c_d)
    E = cubeExpMat(E)
    G = cubeGenMat(G)
    [c, G, Grest, expMat, id] = initialOrderReduction(c, G, Grest, expMat, id, id_, nrGen)
    coeff = leastSquarePolyFunc(x, y, n)
    int = minMaxCubFun(a, b, c, d, l, u)
    int = minMaxDiffCub(a, b, c, d, l, u, f, df)
    int = minMaxDiffQuad(a, b, c, l, u, f, df)
    int = minMaxQuadFun(a, b, c, l, u)
    [c, G, Grest] = quadApproxPolyZono(c, G, Grest, c_a, c_b, c_c)
    [c, G, Grest, expMat, id, d] = reducePolyZono(c, G, Grest, expMat, id, d, nrGen)
    E = squareExpMat(E)
    G_ = squareGenMat(G)
    % --- new ---
    [c, G, Grest] = calcSquared(c1, G1, Grest1, E1, c2, G2, Grest2, E2, isEqual)
    [G_quad] = calcSquaredG(G1, G2, isEqual)
    [G1_ind, G2_ind, G1_ind2, G2_ind2] = calcSquaredGInd(G1, G2, isEqual)
    [E_quad] = calcSquaredE(E1, E2, isEqual)
    der1 = getDerInterval(coeffs, l, u)
    L = minMaxDiffOrder(order, coeffs, l, u, f, der1)
    L = minMaxDiffPoly(p1, p2, l, u)
    coeffs = calcAlternatingDerCoeffs(l, u, order, f, dfs)
    [G_start, G_end, G_ext_start, G_ext_end] = getOrderIndicesG(G, order)
    [Grest_start, Grest_end, Grest_ext_start, Grest_ext_end] = getOrderIndicesGrest(Grest, G, order)
    [c, G, Grest, expMat, id] = restructure(c, G, Grest, expMat, id, id_, nrGen)
    coeff = getChebyshevCoeffs(l, u, order)
    coeff = calcTaylorTanh(order)
    coeff = leastSquareRidgePolyFunc(x, y, n, lambda)
    coeffs = getCoeffsByApproxType(approx_type, order, l, u, f, dfs)
end
end

%------------- END OF CODE --------------
