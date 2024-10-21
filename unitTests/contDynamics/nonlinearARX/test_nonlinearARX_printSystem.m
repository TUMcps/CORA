function res = test_nonlinearARX_printSystem
% test_nonlinearARX_printSystem - unit test function of printSystem
%
% Syntax:
%    res = test_nonlinearARX_printSystem
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
f = @(y,u) [y(1,1) + u(1,1) - y(2,1); ...
              y(3,1) + u(2,1)*cos(y(1,1)); ...
              y(5,1) + u(4,1)*sin(y(1,1))];
dt = 0.25;
sys = nonlinearARX(f,dt,3,2,2);

printSystem(sys)
printSystem(sys,'high')
printSystem(sys,'high',true)
printSystem(sys,'high',false)

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
