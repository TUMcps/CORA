function intVal = transpose(intVal)
% transpose - Overloaded '.'' operator for single operand
%
% Syntax:  
%    intVal = transpose(intVal)
%
% Inputs:
%    intVal - interval object
%
% Outputs:
%    intVal - interval object
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: mtimes

% Author:       Dmitry Grebenyuk
% Written:      07-February-2016
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%infimum
intVal.inf = intVal.inf.';

%supremum
intVal.sup = intVal.sup.';

%------------- END OF CODE --------------
