function res = sum(intVal)
% sum - Overloaded 'sum()' operator for intervals
%
% Syntax:  
%    res = sum(intVal)
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
% Written:      05-August-2016
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% init
res = interval();

%infimum
res.inf = sum(intVal.inf);

%supremum
res.sup = sum(intVal.sup);

%------------- END OF CODE --------------