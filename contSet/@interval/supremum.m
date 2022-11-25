function res = supremum(obj)
% supremum - returns the supremum of an interval
%
% Syntax:  
%    res = sup(obj)
%
% Inputs:
%    obj - interval object
%
% Outputs:
%    res - numerical value
%
% Example: 
%    a = interval([-1 1], [1 2]);
%    b = sup(a)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Matthias Althoff
% Written:      25-June-2015
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

res = obj.sup;

%------------- END OF CODE --------------