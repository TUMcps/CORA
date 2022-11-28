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
%    zB = zonoBundle.generateRandom('Dimension',2);
%    cPZ = conPolyZono(zB)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonoBundle, conZonotope

% Author:       Niklas Kochdumper
% Written:      21-January-2020
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

cPZ = conPolyZono(conZonotope(zB));
    
%------------- END OF CODE --------------