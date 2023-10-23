function res = test_nonlinearSys_display
% test_nonlinearSys_display - unit test for display function
%
% Syntax:
%    res = test_nonlinearSys_display
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

% one-dimensional, no inputs, no outputs
f = @(x,u) x(1)^2;
sys = nonlinearSys(f)

% one-dimensional,  with inputs, no outputs
f = @(x,u) x(1)^2*u(1) - exp(x(1));
sys = nonlinearSys(f)

% one-dimensional, with inputs, linear output
f = @(x,u) x(1)^2*u(1) - exp(x(1));
g = @(x,u) x(1) - u(1);
sys = nonlinearSys(f,g)

% one-dimensional, with inputs, nonlinear output
f = @(x,u) x(1)^2*u(1) - exp(x(1));
g = @(x,u) sin(x(1)) - u(1);
sys = nonlinearSys(f,g)


% multi-dimensional, no inputs
f = @(x,u) [x(1)^2 + x(2);
            x(2) - exp(x(1))];
sys = nonlinearSys(f)

% multi-dimensional, with inputs
f = @(x,u) [x(1)^2 + x(2) - u(1);
            x(2) - exp(x(1)) + u(2)];
sys = nonlinearSys(f)

% multi-dimensional, with inputs, linear output
f = @(x,u) [x(1)^2 + x(2) - u(1);
            x(2) - exp(x(1)) + u(2)];
g = @(x,u) x(2) - x(1);
sys = nonlinearSys(f,g)

% multi-dimensional, with inputs, nonlinear output
f = @(x,u) [x(1)^2 + x(2) - u(1);
            x(2) - exp(x(1)) + u(2)];
g = @(x,u) x(2)*x(1);
sys = nonlinearSys(f,g)


% code executed successfully
res = true;

% ------------------------------ END OF CODE ------------------------------
