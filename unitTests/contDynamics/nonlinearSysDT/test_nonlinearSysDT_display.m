function res = test_nonlinearSysDT_display
% test_nonlinearSysDT_display - unit test for display function
%
% Syntax:
%    res = test_nonlinearSysDT_display
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

dt = 1;

% one-dimensional, no inputs, no outputs
f = @(x,u) x(1)^2;
sys = nonlinearSysDT(f,dt)

% one-dimensional,  with inputs, no outputs
f = @(x,u) x(1)^2*u(1) - exp(x(1));
sys = nonlinearSysDT(f,dt)

% one-dimensional, with inputs, linear output
f = @(x,u) x(1)^2*u(1) - exp(x(1));
g = @(x,u) x(1) - u(1);
sys = nonlinearSysDT(f,dt,g)

% one-dimensional, with inputs, nonlinear output
f = @(x,u) x(1)^2*u(1) - exp(x(1));
g = @(x,u) sin(x(1)) - u(1);
sys = nonlinearSysDT(f,dt,g)


% multi-dimensional, no inputs
f = @(x,u) [x(1)^2 + x(2);
            x(2) - exp(x(1))];
sys = nonlinearSysDT(f,dt)

% multi-dimensional, with inputs
f = @(x,u) [x(1)^2 + x(2) - u(1);
            x(2) - exp(x(1)) + u(2)];
sys = nonlinearSysDT(f,dt)

% multi-dimensional, with inputs, linear output
f = @(x,u) [x(1)^2 + x(2) - u(1);
            x(2) - exp(x(1)) + u(2)];
g = @(x,u) x(2) - x(1);
sys = nonlinearSysDT(f,dt,g)

% multi-dimensional, with inputs, nonlinear output
f = @(x,u) [x(1)^2 + x(2) - u(1);
            x(2) - exp(x(1)) + u(2)];
g = @(x,u) [x(2)*x(1); u(1)*sqrt(x(1))];
sys = nonlinearSysDT(f,dt,g)


% code executed successfully
res = true;

% ------------------------------ END OF CODE ------------------------------
