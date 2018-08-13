function obj = project(obj,dims)
% project - returns a pplPolytope that is projected n the dimesnions
% specified in dims
%
% Syntax:  
%    obj = project(obj,dims)
%
% Inputs:
%    obj - pplPolytope object
%    dims - vector of dimesnions to be projected
%
% Outputs:
%   obj - pplPolytope object
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
% Written:      20-October-2010
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%project dimensions
[C_new d_new] = remove_dim(-obj.C,obj.d,dims);
obj.C = -C_new;
obj.d = d_new;


%------------- END OF CODE --------------