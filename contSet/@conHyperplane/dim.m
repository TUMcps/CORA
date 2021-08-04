function d = dim(hyp)
% isequal - returns the dimension of the given conHyperplane
%
% Syntax:  
%    res = dim(hyp)
%
% Inputs:
%    hyp - conHyperplane object
%
% Outputs:
%    d -  dimension
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
d = max(length(hyp.h.c),size(hyp.C,2));
%------------- END OF CODE --------------