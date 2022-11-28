function P = mptPolytope(hyp)
% mptPolytope - Converts a constrained hyperplane to a polytope
%
% Syntax:  
%    P = mptPolytope(hyp)
%
% Inputs:
%    hyp - conHyperplane object
%
% Outputs:
%    P - mptPolytope object
%
% Example:
%    hyp = conHyperplane(halfspace([1;1],0),[1 0;-1 0],[2;2]);
%    P = mptPolytope(hyp);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: halfspace/mptPolytope

% Author:       Niklas Kochdumper
% Written:      26-Nov-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% conversion
if isempty(hyp.C)
    A = [hyp.h.c';-hyp.h.c'];
    b = [hyp.h.d;-hyp.h.d];       
else
    A = [hyp.h.c';-hyp.h.c';hyp.C];
    b = [hyp.h.d;-hyp.h.d;hyp.d];
end

% instantiate polytope
P = mptPolytope(A,b);

%------------- END OF CODE --------------