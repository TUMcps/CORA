function res = sinh(intVal)
% sinh - Overloaded 'sinh()' operator for intervals
%
% Syntax:  
%    res = sinh(intVal)
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

% [sinh(x_), sinh(x--)].

%------------- BEGIN CODE --------------

res = interval();

res.inf = sinh(intVal.inf);
res.sup = sinh(intVal.sup);

% matrix case
% rand(100, 100)
%time_CORA =
%    5.578342515440901e-004
%time_INTLAB =
%   0.008773613318865


%------------- END OF CODE --------------