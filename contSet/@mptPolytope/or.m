function obj = or(obj,P)
% or - computes union of two mptPolytopes
%
% Syntax:  
%    obj = or(obj,P)
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
% Written:      12-August-2016
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%compute intersection
obj.P = [obj.P P.P];


%------------- END OF CODE --------------