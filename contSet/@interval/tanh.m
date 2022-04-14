function res = tanh(intVal)
% tanh - Overloaded 'tanh()' operator for intervals
%
% x_ is x infimum, x-- is x supremum
%
% [tanh(x_), tanh(x--)].
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
% Last update:  22-February-2016 (DG, the matrix case is rewritten)
% Last revision:---

%------------- BEGIN CODE --------------

res = interval();

res.inf = tanh(intVal.inf);
res.sup = tanh(intVal.sup);

%------------- END OF CODE --------------