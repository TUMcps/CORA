function res = unite_ints(input)
% inites intervals into one
%
% Syntax:  
%    res = unite_ints(input)
%
% Inputs:
%  
%   input - input intervals
%
% Outputs:
%    res - an interval
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: mtimes

% Author:       Dmitry Grebenyuk
% Written:      06-July-2016
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

min_int = infimum(input);
max_int = supremum(input);

res = interval(min(min_int), max(max_int));

%------------- END OF CODE --------------

end
