function d = dim(zB)
% dim - return dimension of zonotope
%
% Syntax:  
%    d = dim(zB)
%
% Inputs:
%    zB - zonoBundle object
%
% Outputs:
%    d - dimension of zB
%
% Example: 
%    zB = zonoBundle.generateRandom(2);
%    d = dim(zB)
%
% Other m-files required: center.m
% Subfunctions: none
% MAT-files required: none
%
% See also: rank.m
%
% Author:        Niklas Kochdumper
% Written:       23-November-2020
% Last update:   ---
% Last revision: ---

%------------- BEGIN CODE --------------

    d = length(center(zB.Z{1}));

%------------- END OF CODE --------------