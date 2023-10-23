function example_manual_nonlinearSys()
% example_manual_nonlinearSys - example from the manual demonstrating the 
% nonlinearSys constructor as defined in the manual
%
% Syntax:
%   example_manual_nonlinearSys()
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

% differential equation f(x,u)
f = @(x,u) [x(2) + u;(1-x(1)^2)*x(2)-x(1)];

% nonlinear system
sys = nonlinearSys(f);

% ------------------------------ END OF CODE ------------------------------
