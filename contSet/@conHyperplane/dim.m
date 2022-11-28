function n = dim(hyp)
% dim - returns the dimension of the ambient space of a constrained
%    hyperplane
%
% Syntax:  
%    n = dim(hyp)
%
% Inputs:
%    hyp - conHyperplane object
%
% Outputs:
%    n - dimension of the ambient space
%
% Example: 
%    hyp = conHyperplane(halfspace([1;1],0),[1 0;-1 0],[2;2]);
%    n = dim(hyp)
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

n = max(length(hyp.h.c),size(hyp.C,2));

%------------- END OF CODE --------------