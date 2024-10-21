function res = test_nonlinearSysDT_printSystem
% test_nonlinearSysDT_printSystem - unit test function of printSystem
%
% Syntax:
%    res = test_nonlinearSysDT_printSystem
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

% Authors:       Tobias Ladner
% Written:       10-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% test print of a simple system
f = @(x,u) [x(1) + u(1);x(2) + u(2)*cos(x(1));x(3) + u(2)*sin(x(1))];
dt = 0.25;
sys = nonlinearSysDT(f,dt);

printSystem(sys)
printSystem(sys,'high')
printSystem(sys,'high',true)
printSystem(sys,'high',false)

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
