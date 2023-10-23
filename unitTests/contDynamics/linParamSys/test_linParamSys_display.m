function res = test_linParamSys_display
% test_linParamSys_display - test the display function
%
% Syntax:
%    res = test_linParamSys_display
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       22-November-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% init system
Ac = [-2 0; 1.5 -3];
Aw = [0 0; 0.5 0];
A_int = intervalMatrix(Ac,Aw);
A_zon = matZonotope(A_int);
B = [1; 1];

% A ... interval matrix, B ... matrix, time-varying parameters
sys = linParamSys(A_int,B,'varParam')

% A ... interval matrix, B ... matrix, time-constant parameters
sys = linParamSys(A_int,B,'constParam')

% A ... matrix zonotope, B ... matrix, time-varying parameters
sys = linParamSys(A_zon,B,'varParam')

% A ... matrix zonotope, B ... matrix, time-constant parameters
sys = linParamSys(A_zon,B,'constParam')

% all checks successful
res = true;

% ------------------------------ END OF CODE ------------------------------
