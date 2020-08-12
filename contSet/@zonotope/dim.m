function d = dim(Z)
% dim - return dimension of zonotope
%
% Syntax:  
%    d = dim(Z)
%
% Inputs:
%    Z - zonotope object
%
% Outputs:
%    d - dimension of Z
%
% Example: 
%    Z = zonotope([zeros(3,1),rand(3,5)]);
%    d = dim(Z)
%
% Other m-files required: center.m
% Subfunctions: none
% MAT-files required: none
%
% See also: rank.m
%
% Author:        Mark Wetzlinger
% Written:       15-Sep-2019
% Last update:   ---
% Last revision: ---

%------------- BEGIN CODE --------------

d = length(center(Z));

%------------- END OF CODE --------------