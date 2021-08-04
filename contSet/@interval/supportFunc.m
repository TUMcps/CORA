function [val,x] = supportFunc(obj,dir,varargin)
% supportFunc - Calculate the upper or lower bound of a interval object
%               along a certain direction
%
% Syntax:  
%    val = supportFunc(obj,dir)
%    val = supportFunc(obj,dir,type)
%
% Inputs:
%    obj - interval object
%    dir - direction for which the bounds are calculated (vector)
%    type - upper or lower bound ('lower' or 'upper')
%
% Outputs:
%    val - bound of the constraind zonotope in the specified direction
%    x - support vector
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: conZonotope/supportFunc

% Author:       Niklas Kochdumper
% Written:      19-November-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

if size(dir,1) == 1 && size(dir,2) >= 1
    % transpose dir
    dir = dir';
end
[val,x] = supportFunc(zonotope(obj),dir,varargin{:});

%------------- END OF CODE --------------