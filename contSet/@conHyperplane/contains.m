function res = contains(hyp,S)
% contains - determines if a constrained hyperplane contains a set or a
%    point
%
% Syntax:  
%    res = contains(hyp,S)
%
% Inputs:
%    hyp - conHyperplane object
%    S - contSet object or single point
%
% Outputs:
%    res - true/false
%
% Example: 
%    hyp = conHyperplane(halfspace([1;1],0),[1 0;-1 0],[2;2]);
%    point = [0;0];
%    contains(hyp,point)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: conHyperplane/isIntersecting

% Author:       Victor Gassmann
% Written:      19-July-2021
% Last update:  25-November-2022 (MW, rename 'contains')
% Last revision:---

%------------- BEGIN CODE --------------

% containment check (input arguments are checked there)
res = contains(mptPolytope(hyp),S);

%------------- END OF CODE --------------