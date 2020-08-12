function d = dim(pZ)
% dim - return dimension of polyZonotope object
%
% Syntax:  
%    d = dim(pZ)
%
% Inputs:
%    pZ - polyZonotope object
%
% Outputs:
%    d - dimension of pZ
%
% Example: 
%    pZ = polyZonotope.generateRandom(4);
%    d = dim(pZ)
%
% Other m-files required: center.m
% Subfunctions: none
% MAT-files required: none
%
% See also: rank.m
%
% Author:        Niklas Kochdumper
% Written:       26-December-2019
% Last update:   ---
% Last revision: ---

%------------- BEGIN CODE --------------

    d = length(pZ.c);

%------------- END OF CODE --------------