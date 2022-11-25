function res = dim(h)
% dim - returns dimension of halfspace
%
% Syntax:  
%    res = dim(h)
%
% Inputs:
%    h - halfspace object
%
% Outputs:
%    res - dimension of halfspace
%
% Example: 
%    h = halfspace([1;1],3);
%    dim(h)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Mark Wetzlinger
% Written:      27-Sep-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

res = length(h.c);

%------------- END OF CODE --------------