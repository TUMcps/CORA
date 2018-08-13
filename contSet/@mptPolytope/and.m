function obj = and(obj,P)
% and - computes intersection of two mptPolytopes
%
% Syntax:  
%    obj = and(obj,P)
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
% Written:      01-February-2011
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%compute intersection
obj.P = obj.P & P.P;


%------------- END OF CODE --------------