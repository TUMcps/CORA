function obj = and(obj,P)
% and - computes intersection of two pplPolytopes
%
% Syntax:  
%    obj = and(obj,P)
%
% Inputs:
%    obj - pplPolytope object
%    P - pplPolytope object
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
% Written:      19-October-2010
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

obj.C = [obj.C; P.C];
obj.d = [obj.d; P.d];

% %compute intersection
% [C,d] = Intersection(-obj.C, -P.C, obj.d, P.d);
% 
% %account for different representation
% obj.C=-C;
% obj.d=d;

%------------- END OF CODE --------------