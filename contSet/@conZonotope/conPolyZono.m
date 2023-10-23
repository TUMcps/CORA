function cPZ = conPolyZono(cZ)
% conPolyZono - convert a constrained zonotope to a constrained polynomial 
%    zonotope
%
% Syntax:
%    cPZ = conPolyZono(cZ)
%
% Inputs:
%    cZ - conZonotope object
%
% Outputs:
%    cPZ - conPolyZono object
%
% Example: 
%    Z = [0 3 0 1;0 0 2 1];
%    A = [1 0 1]; b = 1;
%    cZ = conZonotope(Z,A,b);
%    cPZ = conPolyZono(cZ)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: conZonotope, polyZonotope

% Authors:       Niklas Kochdumper
% Written:       21-January-2020
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

p = size(cZ.G,2);

if ~isempty(cZ.A)
    cPZ = conPolyZono(cZ.c,cZ.G,eye(p),cZ.A,cZ.b,eye(p));
else
    cPZ = conPolyZono(cZ.c,cZ.G,eye(p));
end
    
% ------------------------------ END OF CODE ------------------------------
