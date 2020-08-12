function res = zonoBundle(obj)
% conZonotope - convert an interval into a zonotope bundle object
%
% Syntax:  
%    res = zonoBundle(obj)
%
% Inputs:
%    obj - interval object
%
% Outputs:
%    res - zonoBundle object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: mptPolytope/zonoBundle

% Author:       Niklas Kochdumper
% Written:      26-November-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
    
obj = zonotope(obj);
res = zonoBundle({obj});

%------------- END OF CODE --------------