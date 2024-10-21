function res = test_nonlinearSys_printSystem
% test_nonlinearSys_printSystem - unit test function of printSystem
%
% Syntax:
%    res = test_nonlinearSys_printSystem
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
fun = @(x,u) [x(2); (1-x(1)^2)*x(2)-x(1)];
sys = nonlinearSys('vanDerPol',fun);

printSystem(sys)
printSystem(sys,'high')
printSystem(sys,'high',true)
printSystem(sys,'high',false)

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
