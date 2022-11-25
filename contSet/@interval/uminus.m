function intVal = uminus(intVal)
% uminus - Overloaded '-' operator for single operand
%
% Syntax:  
%    intVal = uplus(intVal)
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

% Author:       Matthias Althoff
% Written:      25-June-2015
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%store values
inf = intVal.inf;
sup = intVal.sup;

%infimum
intVal.inf = -sup;

%supremum
intVal.sup = -inf;

%------------- END OF CODE --------------