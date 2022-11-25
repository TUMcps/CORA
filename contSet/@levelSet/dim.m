function d = dim(obj)
% dim - Returns the dimension of a level set
%
% Syntax:  
%    d = dim(obj)
%
% Inputs:
%    obj - levelSet object
%
% Outputs:
%    d - dimension of the level set
%
% Example: 
%   syms x y
%   eq = x^2 + y^2 - 4;
%   ls = levelSet(eq,[x;y],'==');
%
%   dim(ls)
%
% Other m-files required: center.m
% Subfunctions: none
% MAT-files required: none
%
% See also: zonootpe/dim, interval/dim, levelSet

% Author:       Niklas Kochdumper
% Written:      10-December-2021 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    d = obj.dim;

%------------- END OF CODE --------------