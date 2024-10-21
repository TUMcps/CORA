function res = test_nonlinDASys_printSystem
% test_nonlinDASys_printSystem - unit test function of printSystem
%
% Syntax:
%    res = test_nonlinDASys_printSystem
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
f = @(x,y,u) x(1)+1+u(1);
g = @(x,y,u) (x(1)+1)*y(1) + 2;
sys = nonlinDASys(f,g);

printSystem(sys)
printSystem(sys,'high')
printSystem(sys,'high',true)
printSystem(sys,'high',false)

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
