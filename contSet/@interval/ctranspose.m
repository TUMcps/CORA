function intVal = ctranspose(intVal)
% ctranspose - Overloaded ''' operator for single operand
%
% Syntax:  
%    intVal = ctranspose(intVal)
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
% Written:      14-August-2016
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%infimum
intVal.inf = intVal.inf.';

%supremum
intVal.sup = intVal.sup.';

%------------- END OF CODE --------------
end
