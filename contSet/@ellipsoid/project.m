function E = project(E_in,dims)
% project - Returns an ellipsoid which is projected onto the specified
% dimensions
%
% Syntax:  
%    E = project(E,dims)
%
% Inputs:
%    E - zonotope object
%    dims - projected dimensions
%
% Outputs:
%    E - ellipsoid
%
% Example: 
%    E = ellipsoid([1,0;0,1]);
%    E = project(E, [1 3])
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Victor Gassmann
% Written:      13-March-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
I = eye(length(E_in.Q));
P = I(:,dims);
E = P'*E_in;
%------------- END OF CODE --------------