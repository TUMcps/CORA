function res = isempty(hyp)
% isempty - checks if a constrained hyperplane is the empty set
%
% Syntax:  
%    res = isempty(hyp)
%
% Inputs:
%    hyp - conHyperplane object
%
% Outputs:
%    res - true/false
%
% Example: 
%    hyp1 = conHyperplane(halfspace([1;1],0),[1 0;-1 0],[2;2])
%    isempty(hyp1)
%
%    hyp2 = conHyperplane();
%    isempty(hyp2)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Mark Wetzlinger
% Written:      16-Sep-2019
% Last update:  02-May-2020 (add hyp.d == 0)
% Last revision:---

%------------- BEGIN CODE --------------

res = isempty(hyp.h) && isempty(hyp.C) && (isempty(hyp.d) || hyp.d == 0);

%------------- END OF CODE --------------