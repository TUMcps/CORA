function cPZ = conPolyZono(Z)
% conPolyZono - converts a zonotope to a constrained polynomial zonotope
%
% Syntax:
%    cPZ = conPolyZono(Z)
%
% Inputs:
%    Z - zonotope object
%
% Outputs:
%    cPZ - conPolyZono object
%
% Example: 
%    Z = zonotope([-1;1],[1 3 2 4; 3 2 0 1]);
%    cPZ = conPolyZono(Z)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope, polyZonotope

% Authors:       Niklas Kochdumper
% Written:       21-January-2020
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

cPZ = conPolyZono(polyZonotope(Z));
    
% ------------------------------ END OF CODE ------------------------------
