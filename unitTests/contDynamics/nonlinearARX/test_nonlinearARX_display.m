function res = test_nonlinearARX_display
% test_nonlinearARX_display - unit test for display function
%
% Syntax:
%    res = test_nonlinearARX_display
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

% Authors:       Laura Luetzow
% Written:       15-May-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

dt = 0.1;

% one-dimensional, no inputs, past value
f = @(y,u) 0.5*cos(y(1,1));
sys = nonlinearARX(f,dt, 1, 0, 1)

% one-dimensional, with input, past value
f = @(y,u) 0.5*cos(y(1,1)) + u(1,1) + u(2,1)^2;
sys = nonlinearARX(f,dt, 1, 1, 1)

% one-dimensional, with input, past values
f = @(y,u) 0.5*cos(y(1,1)) + u(1,1) + y(2,1)^2 - sin(u(3,1));
sys = nonlinearARX(f,dt, 1, 1, 2)

% multi-dimensional, no inputs, past value
f = @(y,u) [0.5*cos(y(1,1)); y(2,1)^2];
sys = nonlinearARX(f,dt, 2, 0, 1)

% multi-dimensional, with input, past value
f = @(y,u) [0.5*cos(y(2,1)) + u(1,1) - sin(u(3,1)); y(2,1)-u(2,1)^2];
sys = nonlinearARX(f,dt, 2, 2, 1)

% multi-dimensional, with input, past values
f = @(y,u) [0.5*cos(y(1,1)) + u(1,1) - sin(y(3,1)); y(2,1)^2 - sin(u(3,1) + u(5,1))];
sys = nonlinearARX(f,dt, 2, 2, 2)


% code executed successfully
res = true;

% ------------------------------ END OF CODE ------------------------------
