function res = isempty(hyp)
% isempty - checks if conHyperplane is empty or not
%
% Syntax:  
%    res = isempty(h)
%
% Inputs:
%    hyp - conHyperplane object
%
% Outputs:
%    res - boolean whether hyp is empty or not
%
% Example: 
%    ---
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