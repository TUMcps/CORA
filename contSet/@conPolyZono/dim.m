function d = dim(cPZ)
% dim - return dimension of constrained polynomial zonotope
%
% Syntax:  
%    d = dim(cPZ)
%
% Inputs:
%    cPZ - conPolyZono object
%
% Outputs:
%    d - dimension of cPZ
%
% Example: 
%    cPZ = conPolyZonotope.generateRandom(4);
%    d = dim(cPZ)
%
% Other m-files required: center.m
% Subfunctions: none
% MAT-files required: none
%
% See also: polyZonotope/dim
%
% Author:        Niklas Kochdumper
% Written:       21-January-2021
% Last update:   ---
% Last revision: ---

%------------- BEGIN CODE --------------

    d = length(cPZ.c);

%------------- END OF CODE --------------