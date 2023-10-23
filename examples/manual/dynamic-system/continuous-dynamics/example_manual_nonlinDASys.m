function example_manual_nonlinDASys()
% example_manual_nonlinDASys - example from the manual demonstrating the 
% nonlinDASys constructor as defined in the manual
%
% Syntax:
%   example_manual_nonlinDASys()
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

% differential equation f(x,y,u)
f = @(x,y,u) x + 1 + u;

% constraint equation g(x,y,u)
g = @(x,y,u) (x+1)*y + 2;

% nonlinear differential-algebraic system
sys = nonlinDASys(f,g);

% ------------------------------ END OF CODE ------------------------------
