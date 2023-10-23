function example_manual_nonlinParamSys()
% example_manual_nonlinParamSys - example from the manual demonstrating the 
% nonlinParamSys constructor as defined in the manual
%
% Syntax:
%   example_manual_nonlinParamSys()
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

% differential equation f(x,u,p)
f = @(x,u,p) [x(2) + u;p*(1-x(1)^2)*x(2)-x(1)];

% nonlinear parametric system
sys = nonlinParamSys(f);

% ------------------------------ END OF CODE ------------------------------
