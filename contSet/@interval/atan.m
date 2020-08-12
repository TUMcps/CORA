function res = atan(intVal)
% atan - Overloaded 'atan()' operator for intervals
%
% x_ is x infimum, x-- is x supremum
%
% [atan(x_), atan(x--)].
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
% Last update:  21-February-2016 (DG, the matrix case is rewritten)
% Last revision:---

%------------- BEGIN CODE --------------

res = interval();

res.inf = atan(intVal.inf);
res.sup = atan(intVal.sup);

%------------- END OF CODE --------------