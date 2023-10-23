function cPZ = conPolyZono(pZ)
% conPolyZono - converts a polynomial zonotope to a constrained polynomial
%    zonotope (exactly)
%
% Syntax:
%    cPZ = conPolyZono(pZ)
%
% Inputs:
%    pZ - polyZonotope object
%
% Outputs:
%    cPZ - conPolyZono object
%
% Example: 
%    c = [0;0];
%    G = [2 0 1;0 2 1];
%    GI = [0;0.5];
%    E = [1 0 3;0 1 1];
%    pZ = polyZonotope(c,G,GI,E);
%    cPZ = conPolyZono(pZ)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: polyZonotope

% Authors:       Niklas Kochdumper
% Written:       21-January-2020
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

cPZ = conPolyZono(pZ.c,pZ.G,pZ.E,pZ.GI,pZ.id);
    
% ------------------------------ END OF CODE ------------------------------
