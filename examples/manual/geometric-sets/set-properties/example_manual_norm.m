function example_manual_norm()
% example_manual_norm - example from the manual demonstrating the 
% norm operation as defined in the manual
%
% Syntax:
%   example_manual_norm()
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

% set S
S = zonotope([-0.5 1.5 0; -0.5 0 1.5]);

% norm of the set
res = norm(S,2)

% plot --------------------------------------------------------------------

figure; hold on;
useCORAcolors("CORA:manual")
plot(S)

enlargeAxis(1.2)
title('$\mathcal{S}$','Interpreter','latex');
xlabel('$x_{(1)}$','Interpreter','latex')
ylabel('$x_{(2)}$','Interpreter','latex')

% ------------------------------ END OF CODE ------------------------------
