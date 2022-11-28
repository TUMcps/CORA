function n = dim(pZ)
% dim - returns the dimension of the ambient space of a polynomial zonotope
%
% Syntax:  
%    n = dim(pZ)
%
% Inputs:
%    pZ - polyZonotope object
%
% Outputs:
%    n - dimension of the ambient space
%
% Example: 
%    pZ = polyZonotope.generateRandom('Dimension',4);
%    n = dim(pZ)
%
% Other m-files required: center.m
% Subfunctions: none
% MAT-files required: none
%
% See also: rank.m

% Author:        Niklas Kochdumper
% Written:       26-December-2019
% Last update:   ---
% Last revision: ---

%------------- BEGIN CODE --------------

n = length(pZ.c);

%------------- END OF CODE --------------