function bool_val = iscontained(obj,P)
% iscontained - returns if a polytope P is contained in the polytope
%
% Syntax:  
%    obj = iscontained(obj,P)
%
% Inputs:
%    obj - pplPolytope object
%    P - pplPolytope object
%
% Outputs:
%   bool_val - bool value
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
% Written:      18-November-2010
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------


%check for containment
bool_val = is_contained(-obj.C, -P.C, obj.d, P.d);


%------------- END OF CODE --------------