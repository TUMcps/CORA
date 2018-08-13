function obj = enlarge(obj,factor)
% enlarge - Enlarges an interval object around its center
%
% Syntax:  
%    obj = enlarge(obj,factor)
%
% Inputs:
%    obj - interval object
%    factor - enlarging factor (scalar or column vector)
%
% Outputs:
%    obj - enlarged interval object
%
% Example: 
%    I = interval([1 2; -1 1]);
%    I = enlatge(I,2);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also:

% Author:       Matthias Althoff
% Written:      22-July-2016 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%get center and radius
c = 0.5*(obj.inf + obj.sup);
r = 0.5*(obj.sup - obj.inf);

%enlarged intervals
obj.inf = c-r.*factor; 
obj.sup = c+r.*factor;

%------------- END OF CODE --------------