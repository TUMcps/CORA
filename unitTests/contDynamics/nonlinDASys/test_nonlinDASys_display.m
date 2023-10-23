function res = test_nonlinDASys_display
% test_nonlinDASys_display - unit test for display function
%
% Syntax:
%    res = test_nonlinDASys_display
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

% basic
f = @(x,y,u) x(1) - y(1) + u(1);
g = @(x,y,u) y(1) + x(1)*u(1);
sys = nonlinDASys(f,g)

% multi-dimensional
f = @(x,y,u) [x(1) - y(1) + u(1); y(1)^2 - x(2)^2*u(2)];
g = @(x,y,u) y(1) + x(1)*u(1);
h = @(x,y,u) x(1)*y(1);
sys = nonlinDASys(f,g,h)

% no constraints
f = @(x,y,u) [x(1) - u(1); x(1)^2 - x(2)^2*u(2)];
g = @(x,y,u) x(1) + u(1);
h = @(x,y,u) x(1)*u(1);
sys = nonlinDASys(f,g,h)

% multi-dimensional
f = @(x,y,u) [x(1) + y(1) - u(1); log(x(3)^2) - x(2)^2*u(2); 0];
g = @(x,y,u) [x(1)*sin(u(1)); y(2) - x(1)];
h = @(x,y,u) [x(1)*u(1); 1; sqrt(x(3))*cos(y(1))];
sys = nonlinDASys(f,g,h)


% code executed successfully
res = true;

% ------------------------------ END OF CODE ------------------------------
