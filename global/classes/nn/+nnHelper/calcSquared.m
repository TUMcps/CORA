function [c, G, GI] = calcSquared(c1, G1, GI1, E1, c2, G2, GI2, E2, isEqual)
% calcSquared - computes pZ_1 * pZ_2 with pZ_1, pZ_2 being multiplicatives 
%    of a one-dimensional polyZonotope pZ:
%    pZ_1 * pZ_2 =
%         = (c1 + G1 + GI1)*(c2 + G2 + GI2)
%         = c1*c2 + c1*G2 + c1*GI2
%           + G1*c2 + G1*G2 + G2*GI2
%           + GI1*c2 + GI1*G2 + GI1*GI2
%         = (c1*c2 + sum(0.5[GI1*GI2](:, I)))
%           % = c
%           + (c1*G2 + G1*c2 + G1*G2)
%           % = G
%           + (c1*GI2 + GI2*c2 + 0.5[GI1*GI2](:, I) + [GI1*GI2](:, ~I) + G1*GI2 + GI2*G2)
%           % = GI
%           % with I ... indices of generators that need shifting
%
% Note: (Half of) [GI1*GI2](:, I) appears in c & GI:
%   Squaring independent Generators need shifting
%   as only positive part is used afterwards.
%   -> this will be corrected when adding the terms together
%
% Syntax:
%    res = nnHelper.calcSquared(c1, G1, GI1, c2, G2, GI2, i1, i2, Es)
%
% Inputs:
%    c1 - center of first polyZonotope
%    G1 - dependent generators of first polyZonotope
%    GI1 - independent generators of first polyZonotope
%    c2 - center of seconds polyZonotope
%    G2 - dependent generators of seconds polyZonotope
%    GI2 - independent generators of seconds polyZonotope
%    i1 - index of first polyZonotope
%    i2 - index of second polyZonotope
%    Es - cell array holding all exponential matrices
%
% Outputs:
%    [c, G, GI] - of pZ^(i1+i2)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: polyZonotope/quadMap, nnHelper/calcSquaredG, nnHelper/calcSquaredE

% Authors:       Tobias Ladner
% Written:       28-March-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if nargin < 9
    throw(CORAerror('CORA:notEnoughInputArgs',9));
end

G_quad = nnHelper.calcSquaredG(G1, G2, isEqual);
GI_quad = nnHelper.calcSquaredG(GI1, GI2, isEqual);

% construct squared parameters
c = c1 * c2;

% See Note
if isEqual
    r = length(GI1);
    GI_quad(:, 1:r) = 0.5 * GI_quad(:, 1:r);
    c = c + sum(GI_quad(:, 1:r));
end

% the same principle applies to G1, G2:
% if all exponents are even and they become independent generators
% after multiplying them with GI2, GI1, respectively.
% except center does not need shifting as
% G1, G2 only scale the independent generators of GI2, GI1.

% G1 * GI2
even_indices = all(mod(E1, 2) == 0, 1);
G1_ind = G1; % copy by value
G1_ind(:, even_indices) = 0.5 * G1_ind(:, even_indices);
G1GI2 = nnHelper.calcSquaredG(G1_ind, GI2);

% GI1 * G2
even_indices = all(mod(E2, 2) == 0, 1);
G2_ind = G2; % copy by value
G2_ind(:, even_indices) = 0.5 * G2_ind(:, even_indices);
GI1G2 = nnHelper.calcSquaredG(GI1, G2_ind);

G = [G1 * c2, c1 * G2, G_quad];

if isEqual
    GI = [GI1 * c2, c1 * GI1, GI_quad, 2 * GI1G2];
else
    GI = [GI1 * c2, c1 * GI2, GI_quad, G1GI2, GI1G2];
end

end

% ------------------------------ END OF CODE ------------------------------
