function res = sinh(intVal)
% sinh - Overloaded 'sinh()' operator for intervals
%
% x_ is x infimum, x-- is x supremum
%
% [sinh(x_), sinh(x--)].
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
% Last update:  21-February-2016 (DG, the matrix case is rewritten)
% Last revision:---

%------------- BEGIN CODE --------------

res = interval();

res.inf = sinh(intVal.inf);
res.sup = sinh(intVal.sup);


%------------- END OF CODE --------------