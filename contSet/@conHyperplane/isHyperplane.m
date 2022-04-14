function res = isHyperplane(hyp)
% isequal - checks whether the conHyperplane object is representable as a
% simple hyperplane or not
%
% Syntax:  
%    res = isHyperplane(hyp)
%
% Inputs:
%    hyp - conHyperplane object
%
% Outputs:
%    res - boolean indicating whether hyp is a simple hyperplane or not
%
% Example: 
%    ---
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
% if is one-dimensional and consistent (checked in constructor), always
% true
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
    if supportFunc(hyp,B(:,1),'upper')<inf ||...
       supportFunc(hyp,B(:,1),'lower')>-inf
        res = false;
        return;
    end
end
res = true;
%------------- END OF CODE --------------