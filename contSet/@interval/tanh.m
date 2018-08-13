function res = tanh(intVal)
% tanh - Overloaded 'tanh()' operator for intervals
%
% Syntax:  
%    res = tanh(intVal)
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
% See also: mtimes

% Author:       Matthias Althoff
% Written:      05-February-2016
% Last update:  22-February-2016 the matrix case is rewritten  (Dmitry Grebenyuk)
% Last revision:---

% x_ is x infimum, x-- is x supremum

% [tanh(x_), tanh(x--)].

%------------- BEGIN CODE --------------

res = interval();

res.inf = tanh(intVal.inf);
res.sup = tanh(intVal.sup);

% matrix case
% rand(100, 100)
%time_CORA =
%    3.929944670515823e-004
%time_INTLAB =
%   0.005584974001024

%------------- END OF CODE --------------