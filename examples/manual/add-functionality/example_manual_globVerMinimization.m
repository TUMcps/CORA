function example_manual_globVerMinimization()
% example_manual_globVerMinimization - example from the manual demontrating the
% globVerMinimization operation as defined in the manual
%
% Syntax:
%   example_manual_globVerMinimization()
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

% Authors:       Tobias Ladner
% Written:       27-September-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% function f
f = @(x) (1.5 - x(1)*(1-x(2))).^2 + ...
    (2.25 - x(1)*(1-x(2)^2))^2 + ...
    (2.625 - x(1)*(1-x(2)^3))^2;

% domain D
D = interval([-4.5;-4.5],[4.5;4.5]);

% verified global optimization
[val,xOpt,domOpt] = globVerMinimization(f,D,1e-5);

% ------------------------------ END OF CODE ------------------------------
