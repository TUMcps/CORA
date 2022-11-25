function res = in(H,obj)
% isIntersecting - determines if obj is contained in H
%
% Syntax:  
%    res = in(H,obj)
%
% Inputs:
%    H - conHyperplane object
%    obj - contSet object
%
% Example: 
%    hp = conHyperplane(halfspace([1;1],0),[1 0;-1 0],[2;2]);
%    point = [0;0];
%    in(hp,point)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: conHyperplane/isIntersecing

% Author:       Victor Gassmann
% Written:      19-July-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
res = in(mptPolytope(H),obj);
%------------- END OF CODE --------------