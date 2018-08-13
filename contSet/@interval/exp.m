function res = exp(intVal)
% exp - Overloaded 'exp()' operator for intervals
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

% x_ is x infimum, x-- is x supremum

% [exp(x_), exp(x--)].

%------------- BEGIN CODE --------------

%exponential function is monotonic
res = interval(exp(intVal.inf), exp(intVal.sup));

% matrix case
% rand(100, 100)
%time_CORA =
%    6.127413732464556e-004
%time_INTLAB =
%   0.004009680390026

%------------- END OF CODE --------------