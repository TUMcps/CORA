function res = isempty(hs)
% isempty - checks if a halfspace is the empty set
%
% Syntax:  
%    res = isempty(hs)
%
% Inputs:
%    hs - halfspace object
%
% Outputs:
%    res - true/false
%
% Example: 
%    hs = halfspace([1 1],2);
%    isempty(hs)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Mark Wetzlinger
% Written:      16-Sep-2019
% Last update:  02-May-2020 (adapt to new default values in constructor)
% Last revision:---

%------------- BEGIN CODE --------------

res = isempty(hs.c) && (isempty(hs.d) || hs.d == 0);

%------------- END OF CODE --------------