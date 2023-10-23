function E = enlarge(E,factor)
% enlarge - enlarge ellipsoid by a factor
%
% Syntax:
%    E = enlarge(E,factor)
%
% Inputs:
%    E - ellipsoid object
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

% Authors:       Mark Wetzlinger
% Written:       15-September-2019
% Last update:   04-July-2022 (VG, input checks)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check input arguments
inputArgsCheck({{E,'att','ellipsoid','scalar'}; ...
                {factor,'att','numeric',{'scalar','nonnegative'}}});

% enlarge shape matrix
E = ellipsoid(factor^2*E.Q,center(E),E.TOL);

% ------------------------------ END OF CODE ------------------------------
