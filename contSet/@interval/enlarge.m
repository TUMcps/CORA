function I = enlarge(I,factor)
% enlarge - Enlarges an interval object around its center
%
% Syntax:  
%    obj = enlarge(I,factor)
%
% Inputs:
%    I - interval object
%    factor - enlarging factor (scalar or column vector)
%
% Outputs:
%    I - enlarged interval object
%
% Example: 
%    I = interval([-1;-2],[3;4]);
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
c = center(I);
r = rad(I);

%enlarged intervals
I.inf = c-r.*factor; 
I.sup = c+r.*factor;

%------------- END OF CODE --------------