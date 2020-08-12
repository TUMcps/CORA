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
%    I = enlarge(I,2);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also:

% Author:       Matthias Althoff
% Written:      22-July-2016 
% Last update:  28-Aug-2019
% Last revision:---

%------------- BEGIN CODE --------------

%get center and radius
c = center(obj);
r = rad(obj);

%enlarged intervals
obj.inf = c-r.*factor; 
obj.sup = c+r.*factor;

%------------- END OF CODE --------------