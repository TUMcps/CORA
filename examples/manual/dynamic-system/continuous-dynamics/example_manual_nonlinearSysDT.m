function example_manual_nonlinearSysDT()
% example_manual_nonlinearSysDT - example from the manual demonstrating the 
% nonlinearSysDT constructor as defined in the manual
%
% Syntax:
%   example_manual_nonlinearSysDT()
%
% Inputs:
%    -
%
% Outputs:
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also:

% Authors:        Tobias Ladner
% Written:       27-September-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% equation f(x,u)
f = @(x,u) [x(1) + u(1); ...
    x(2) + u(2)*cos(x(1)); ...
    x(3) + u(2)*sin(x(1))];

% sampling time
dt = 0.25;

% nonlinear discrete-time system
sys = nonlinearSysDT(f,dt);

% ------------------------------ END OF CODE ------------------------------
