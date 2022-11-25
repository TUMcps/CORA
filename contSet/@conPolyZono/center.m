function c = center(cPZ)
% center - returns the center of a constrained polynomial zonotope
%
% Syntax:  
%    c = center(cPZ)
%
% Inputs:
%    cPZ - conPolyZonotope object
%
% Outputs:
%    c - center of the constrained polynomial zonotope
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/center, quadZonotope/center

% Author:       Niklas Kochdumper
% Written:      14-August-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    c = cPZ.c;

%------------- END OF CODE --------------