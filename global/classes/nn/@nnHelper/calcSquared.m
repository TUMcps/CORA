function [c, G, Grest] = calcSquared(c1, G1, Grest1, E1, c2, G2, Grest2, E2, isEqual)
% calcSquared - computes pZ^i1 * pZ^i2
%    with pZ^i1, pZ^i2 being multiplicatives of a polyZonotope pZ
%    -> (c1 + G1 + Grest1)*(c2 + G2 + Grest2) =
%         = c1*c2 + c1*G2 + c1*Grest2
%           + G1*c2 + G1*G2 + G2*Grest2
%           + Grest1*c2 + Grest1*G2 + Grest1*Grest2
%         = (c1*c2 + sum(0.5[Grest1*Grest2](:, I)))
%           % = c
%           + (c1*G2 + G1*c2 + G1*G2)
%           % = G
%           + (c1*Grest2 + Grest2*c2 + 0.5[Grest1*Grest2](:, I) + [Grest1*Grest2](:, ~I) + G1*Grest2 + Grest2*G2)
%           % = Grest
%           % with I ... indices of generators that need shifting
%
% Note: (Half of) [Grest1*Grest2](:, I) appears in c & Grest:
%   Squaring independent Generators need shifting
%   as only positive part is used afterwards.
%   -> this will be corrected when adding the terms together
%
% Syntax:
%    res = nnHelper.calcSquared(c1, G1, Grest1, c2, G2, Grest2, i1, i2, Es)
%
% Inputs:
%    c1 - center of first polyZonotope
%    G1 - dependent generators of first polyZonotope
%    Grest1 - independent generators of first polyZonotope
%    c2 - center of seconds polyZonotope
%    G2 - dependent generators of seconds polyZonotope
%    Grest2 - independent generators of seconds polyZonotope
%    i1 - index of first polyZonotope
%    i2 - index of second polyZonotope
%    Es - cell array holding all exponential matrices
%
% Outputs:
%    [c, G, Grest] - of pZ^(i1+i2)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: polyZonotope/quadMap, nnHelper/calcSquaredG, nnHelper/calcSquaredE
%    Generalization of old nnHelper/quadApproxPolyZono & cubApproxPolyZono

% Author:       Tobias Ladner
% Written:      28-March-2022
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

if nargin < 9
    throw(CORAerror('CORA:notEnoughInputArgs',9));
end

G_quad = nnHelper.calcSquaredG(G1, G2, isEqual);
Grest_quad = nnHelper.calcSquaredG(Grest1, Grest2, isEqual);

% construct squared parameters
c = c1 * c2;

% See Note
if isEqual
    r = length(Grest1);
    Grest_quad(:, 1:r) = 0.5 * Grest_quad(:, 1:r);
    c = c + sum(Grest_quad(:, 1:r));
end

% the same principle applies to G1, G2:
% if all exponents are even and they become independent generators
% after multiplying them with Grest2, Grest1, respectively.
% except center does not need shifting as
% G1, G2 only scale the independent generators of Grest2, Grest1.

% G1 * Grest2
even_indices = all(mod(E1, 2) == 0, 1);
G1_ind = G1; % copy by value
G1_ind(:, even_indices) = 0.5 * G1_ind(:, even_indices);
G1Grest2 = nnHelper.calcSquaredG(G1_ind, Grest2);

% Grest1 * G2
even_indices = all(mod(E2, 2) == 0, 1);
G2_ind = G2; % copy by value
G2_ind(:, even_indices) = 0.5 * G2_ind(:, even_indices);
Grest1G2 = nnHelper.calcSquaredG(Grest1, G2_ind);

G = [G1 * c2, c1 * G2, G_quad];

if isEqual
    Grest = [Grest1 * c2, c1 * Grest1, Grest_quad, 2 * Grest1G2];
else
    Grest = [Grest1 * c2, c1 * Grest2, Grest_quad, G1Grest2, Grest1G2];
end

end

%------------- END OF CODE --------------