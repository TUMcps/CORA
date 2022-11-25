function res = project(obj, dim)
% project - project interval onto given dimensions
%
% Syntax:  
%    res = project(obj, dim)
%
% Inputs:
%    obj - interval object
%    dim - dimensions for projection
%
% Outputs:
%    res - interval
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Author:       Mark Wetzlinger
% Written:      16-Sep-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

infi = infimum(obj);
supr = supremum(obj);
res = interval(infi(dim),supr(dim));

%------------- END OF CODE --------------