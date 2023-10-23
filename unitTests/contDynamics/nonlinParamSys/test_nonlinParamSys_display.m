function res = test_nonlinParamSys_display
% test_nonlinParamSys_display - unit test for display function
%
% Syntax:
%    res = test_nonlinParamSys_display
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

% one-dimensional, no inputs, no parameters, no outputs
f = @(x,u,p) x(1)^2;
sys = nonlinParamSys(f)

% one-dimensional, with inputs, no parameters, no outputs
f = @(x,u,p) x(1)^2*u(1) - exp(x(1));
sys = nonlinParamSys(f)

% one-dimensional, with inputs, with parameters, no outputs
f = @(x,u,p) x(1)^2*u(1) - exp(p(1));
sys = nonlinParamSys(f)

% one-dimensional, with inputs, with parameters, linear output
f = @(x,u,p) x(1)^2*u(1) - exp(p(1));
g = @(x,u,p) x(1) - u(1);
sys = nonlinParamSys(f,g)

% one-dimensional, with inputs, with parameters, nonlinear output
f = @(x,u,p) x(1)^2*u(1) - exp(p(1));
g = @(x,u,p) sin(x(1)) - u(1)*p(1);
sys = nonlinParamSys(f,g)


% multi-dimensional, no inputs, no parameters
f = @(x,u,p) [x(1)^2 + x(2);
            x(2) - exp(x(1))];
sys = nonlinParamSys(f)

% multi-dimensional, with inputs, no parameters
f = @(x,u,p) [x(1)^2 + x(2) - u(1);
            x(2) - exp(x(1)) + u(2)];
sys = nonlinParamSys(f)

% multi-dimensional, with inputs, with parameters
f = @(x,u,p) [x(1)^2 + x(2) - u(1);
            x(2)*p(1) - exp(x(1)) + p(2)];
sys = nonlinParamSys(f)

% multi-dimensional, with inputs, no parameters, linear output
f = @(x,u,p) [x(1)^2 + x(2) - u(1);
            x(2) - exp(x(1)) + u(2)];
g = @(x,u,p) x(2) - x(1);
sys = nonlinParamSys(f,g)

% multi-dimensional, with inputs, with parameters, linear output
f = @(x,u,p) [x(1)^2 + x(2) - u(1);
            x(2)*p(1) - exp(x(1)) + u(2)];
g = @(x,u,p) x(2) - p(1);
sys = nonlinParamSys(f,g)

% multi-dimensional, with inputs, no parameters, nonlinear output
f = @(x,u,p) [x(1)^2 + x(2) - u(1);
            x(2) - exp(x(1)) + u(2)];
g = @(x,u,p) x(2)*x(1);
sys = nonlinParamSys(f,g)

% multi-dimensional, with inputs, with parameters, nonlinear output
f = @(x,u,p) [x(1)^2 + x(2) - u(1);
            x(2)*p(1) - exp(x(1)) + p(2)];
g = @(x,u,p) x(2)*p(1);
sys = nonlinParamSys(f,g)


% code executed successfully
res = true;

% ------------------------------ END OF CODE ------------------------------
