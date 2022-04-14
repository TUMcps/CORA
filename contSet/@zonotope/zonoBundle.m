function res = zonoBundle(obj)
% conZonotope - convert a zonotope object into a zonotope bundle object
%
% Syntax:  
%    res = zonoBundle(obj)
%
% Inputs:
%    obj - zonotope object
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
    
    res = zonoBundle({obj});

%------------- END OF CODE --------------