function E = enlarge(E,factor)
% enlarge - enlarge ellipsoid E by factor
%
% Syntax:  
%    E = enlarge(E,factor)
%
% Inputs:
%    E - Ellipsoid object
%    factor - enlargement factor (scalar)
%
% Outputs:
%    E - enlarged E
%
% Example: 
%    E = ellipsoid([1 0; 0 2],[0;0]);
%    factor = 2;
%    E = enlarge(E,factor);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---
%
% Author:       Mark Wetzlinger
% Written:      15-Sep-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

E = ellipsoid(factor^2*E.Q,center(E),E.TOL);

%------------- END OF CODE --------------