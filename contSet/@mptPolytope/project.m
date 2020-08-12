function obj = project(obj,dims)
% project - projects a mptPolytope onto a set of dimensions
%
% Syntax:  
%    res = project(obj,dims)
%
% Inputs:
%    obj - mptPolytope object
%    dims - vector of dimensions
%
% Outputs:
%    res - projected mptPolytope object
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Matthias Althoff
% Written:      14-November-2016
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%compute projection
obj.P = projection(obj.P, dims);

%------------- END OF CODE --------------