function res = isempty(h)
% isempty - checks if halfspace is empty or not
%
% Syntax:  
%    res = isempty(h)
%
% Inputs:
%    h - halfspace object
%
% Outputs:
%    res - boolean whether h is empty or not
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
% Last update:  02-May-2020 (adapt to new default values in constructor)
% Last revision:---

%------------- BEGIN CODE --------------

res = isempty(h.c) && (isempty(h.d) || h.d == 0);

%------------- END OF CODE --------------