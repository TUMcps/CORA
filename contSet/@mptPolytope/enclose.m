function obj = enclose(obj,P)
% enclose - computes the convex hull of two mptPolytopes
%
% Syntax:  
%    obj = enclose(obj,P)
%
% Inputs:
%    obj - mptPolytope object
%    P - mptPolytope object
%
% Outputs:
%   obj - mptPolytope object
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
% Written:      02-February-2011
% Last update:  12-August-2016
% Last revision:---

%------------- BEGIN CODE --------------

try %MPT3
    obj.P = convexHull(obj.P, P.P);
catch %MPT2
    %compute convex hull
    obj.P = hull([obj.P P.P]);
end


%------------- END OF CODE --------------