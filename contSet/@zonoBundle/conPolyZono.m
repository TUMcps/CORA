function cPZ = conPolyZono(zB)
% conPolyZono - convert a zonotope bundle to a constrained polynomial
%    zonotope
%
% Syntax:
%    cPZ = conPolyZono(zB)
%
% Inputs:
%    zB - zonoBundle object
%
% Outputs:
%    cPZ - conPolyZono object
%
% Example: 
%    Z1 = zonotope([1;-1],[2 -1 3; 0 1 -1]);
%    Z2 = Z1 + [1;0];
%    zB = zonoBundle({Z1,Z2});
%    cPZ = conPolyZono(zB)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonoBundle, conZonotope

% Authors:       Niklas Kochdumper
% Written:       21-January-2020
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

cPZ = conPolyZono(conZonotope(zB));
    
% ------------------------------ END OF CODE ------------------------------
