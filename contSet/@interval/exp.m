function res = exp(intVal)
% exp - Overloaded 'exp()' operator for intervals
%
% x_ is x infimum, x-- is x supremum
%
% [exp(x_), exp(x--)].
%
% Syntax:  
%    res = exp(intVal)
%
% Inputs:
%    intVal - interval object
%
% Outputs:
%    res - interval object
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: interval

% Author:       Matthias Althoff
% Written:      25-June-2015
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%exponential function is monotonic
res = interval(exp(intVal.inf), exp(intVal.sup));

%------------- END OF CODE --------------