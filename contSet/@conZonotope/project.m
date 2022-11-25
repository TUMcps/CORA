function obj = project(obj,projDim)
% project - project a constrained zonotope object to a subspace
%
% Syntax:  
%    obj = project(obj,projDim)
%
% Inputs:
%    obj - conZonotope object
%    projDim - dimensions of the projection
%
% Outputs:
%    obj - contrained zonotope object
%
% Example: 
%    Z = [0 1 0 1 0;0 1 2 -1 0;0 0 0 0 1];
%    A = [-2 1 -1 0];
%    b = 2;
%    cZono = conZonotope(Z,A,b);
%    projZono = project(cZono,[1,2]);
%    
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Niklas Kochdumper, Mark Wetzlinger
% Written:      11-May-2018
% Last update:  14-March-2021 (MW, input args handling, empty set)
% Last revision:---

%------------- BEGIN CODE --------------

try
    % try projection first, as isempty() can be quite costly
    obj.Z = obj.Z(projDim,:);
catch ME
    % now check if conZonotope was empty
    if isempty(obj)
        [msg,id] = errEmptySet();
        error(id,msg);
    elseif any(projDim > dim(obj))
        error("Dimensions for projection exceed dimension of conZonotope!");
    else
        rethrow(ME);
    end
end

%------------- END OF CODE --------------