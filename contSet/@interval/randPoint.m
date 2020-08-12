function p = randPoint(obj)
% randPoint - computes random point in interval obj
%
% Syntax:  
%    p = randPoint(obj)
%
% Inputs:
%    obj - interval object
%
% Outputs:
%    p - random point in interval
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Mark Wetzlinger
% Written:      17-Sep-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

c = center(obj);
r = rad(obj);

p = c + (-1 + 2 * rand(length(r),1)) .* r;


%------------- END OF CODE --------------