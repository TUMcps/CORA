function res = atan(intVal)
% atan - Overloaded 'atan()' operator for intervals
%
% Syntax:  
%    res = atan(intVal)
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
% Last update:  21-February-2016 the matrix case is rewritten  (Dmitry Grebenyuk)
% Last revision:---

% x_ is x infimum, x-- is x supremum

% [atan(x_), atan(x--)].

%------------- BEGIN CODE --------------

res = interval();

res.inf = atan(intVal.inf);
res.sup = atan(intVal.sup);

% matrix case
% rand(100, 100)
%time_CORA =
%    3.998627914054368e-004
%time_INTLAB =
%   0.004031745868841

%------------- END OF CODE --------------