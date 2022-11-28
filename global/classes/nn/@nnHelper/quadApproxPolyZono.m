function [c, G, Grest] = quadApproxPolyZono(c, G, Grest, c_a, c_b, c_c)
% quadApproxPolyZono - evaluate a quadratic function c_a*x^2 + c_b*x + c_c 
%    on a polynomial zonotope
%
% Syntax:
%     [c, G, Grest] = quadApproxPolyZono(c, G, Grest, c_a, c_b, c_c)
%
% Inputs:
%    c - center of polyZonotope in a dimension
%    G - corresponding dep. generator of polyZonotope as row vector
%    Grest - corresponding indep. generator of polyZonotope as row vector
%    c_a-c_c - polynomial coefficients
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

q = length(Grest);

% compute (c + G + Grest)*(c + G + Grest) =
%         = c*c + 2*c*G + G*G + 2*c*Grest + 2*G*Grest + Grest*Grest
G_quad = nnHelper.squareGenMat(G);
Grest_quad = nnHelper.squareGenMat(Grest);
Grest_quad(:, 1:q) = 0.5 * Grest_quad(:, 1:q);

% construct resulting polynomial zonotope
Grest = [c_b * Grest + 2 * c_a * c * Grest, c_a * Grest_quad, ...
    2 * c_a * reshape(G'*Grest, 1, [])];
G = [c_b * G + 2 * c_a * c * G, c_a * G_quad];
c = c_a * (c^2 + sum(Grest_quad(:, 1:q))) + c_b * c + c_c;
end

%------------- END OF CODE --------------