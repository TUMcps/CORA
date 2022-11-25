function res = asinh(intVal)
% asinh - Overloaded 'asinh()' operator for intervals
%
% x_ is x infimum, x-- is x supremum
%
% [asinh(x_), asinh(x--)].
%
% Syntax:  
%    res = asinh(intVal)
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
% Written:      12-February-2016
% Last update:  21-February-2016 (DG, the matrix case is rewritten)
% Last revision:---

%------------- BEGIN CODE --------------

res = interval();

res.inf = asinh(intVal.inf);
res.sup = asinh(intVal.sup);

%------------- END OF CODE --------------