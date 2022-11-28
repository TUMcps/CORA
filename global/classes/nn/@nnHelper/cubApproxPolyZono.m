function [c, G, Grest] = cubApproxPolyZono(c, G, Grest, c_a, c_b, c_c, c_d)
% cubApproxPolyZono - evaluate a cubic function c_a*x^3 + c_b*x^2 + c_c*x + c_d on a polynomial
%    zonotope
%
% Syntax:
%    [c, G, Grest] = nnHelper.cubApproxPolyZono(c, G, Grest, c_a, c_b, c_c, c_d)
%
% Inputs:
%    c - center of polyZonotope in a dimension
%    G - corresponding dep. generator of polyZonotope as row vector
%    Grest - corresponding indep. generator of polyZonotope as row vector
%    c_a-c_d - polynomial coefficients
%
% Outputs:
%    [c, G, Grest] - resulting polynomial zonotope
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Author:        Niklas Kochdumper, Tobias Ladner
% Written:       17-September-2021
% Last update:   ---
% Last revision: 28-March-2022 (TL)

%------------- BEGIN CODE --------------

h = length(G);
q = length(Grest);

% compute (c + G + Grest)^3 = (c^3 + 3*c^2*G + 3*c*G*G + G*G*G) + ...
%  3*c^2*Grest + 6*c*G*Grest + 3*c*Grest*Grest + 3*G*Grest*G + ...
%  3*Grest*Grest*G + G_r*G_r*G_r
G_quad = nnHelper.squareGenMat(G);
G_quad_ = G_quad;
G_quad_(:, 1:h) = 0.5 * G_quad_(:, 1:h);
c_quad = 0.5 * sum(G_quad(:, 1:h));
G_cub = nnHelper.cubeGenMat(G);
Grest_quad = nnHelper.squareGenMat(Grest);
Grest_quad(:, 1:q) = 0.5 * Grest_quad(:, 1:q);
crest_quad = 0.5 * sum(Grest_quad(:, 1:q));
Grest_cub = nnHelper.cubeGenMat(Grest);

% construct resulting polynomial zonotope
Grest = [(c_a * 3 * c^2 + c_b * 2 * c + c_c + c_a * 3 * c_quad) * Grest, ...
    (c_a * 6 * c + c_b * 2) * reshape(G'*Grest, 1, []), ...
    (c_a * 3 * c + c_b) * Grest_quad, ...
    (c_a * 3) * reshape(Grest'*G_quad_, 1, []), ...
    (c_a * 3) * reshape(G'*Grest_quad, 1, []), c_a * Grest_cub];
G = [(c_a * 3 * c^2 + c_b * 2 * c + c_a * 3 * crest_quad + c_c) * G, ...
    (c_a * 3 * c + c_b) * G_quad, c_a * G_cub];
c = c_a * c^3 + c_b * c^2 + c_c * c + c_d + (c_a * 3 * c + c_b) * crest_quad;
end

%------------- END OF CODE --------------
