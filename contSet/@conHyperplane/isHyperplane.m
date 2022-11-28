function res = isHyperplane(hyp)
% isHyperplane - checks whether a constrained hyperplane can be represented
%    as a simple hyperplane
%
% Syntax:  
%    res = isHyperplane(hyp)
%
% Inputs:
%    hyp - conHyperplane object
%
% Outputs:
%    res - true/false
%
% Example: 
%    hyp = conHyperplane(halfspace([1;1],0),[1 0;-1 0],[2;2]);
%    isHyperplane(hyp)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Victor Gassmann
% Written:      16-March-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% always true if conHyperplane is 1D and consistent (checked in constructor)
if dim(hyp)==1
    res = true;
    return;
end
if isempty(hyp.C) || (all(all(hyp.C==0)) && all(hyp.d>=0))
    res = true;
    return;
end

% check if C is unbounded for all x on the hyperplane by checking if
% support function is unbounded for all 2*(n-1) directions ('upper' and 'lower')
c = hyp.h.c/norm(hyp.h.c);
n = dim(hyp);
% null space has exactly n-1 vectors
B = null(c');
for i=1:n-1
    if supportFunc(hyp,B(:,1),'upper') < Inf || ...
       supportFunc(hyp,B(:,1),'lower') > -Inf
        res = false;
        return;
    end
end

res = true;

%------------- END OF CODE --------------