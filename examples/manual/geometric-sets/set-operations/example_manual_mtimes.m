function example_manual_mtimes()
% example_manual_mtimes - example from the manual demonstrating the linear
% map of a set as defined in the manual
%
% Syntax:
%   example_manual_mtimes()
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

% set and matrix
S = zonotope([0 1 1 0; 0 1 0 1]);
M = [1 0; -1 0.5];

% linear transformation
res = M*S;

% plot --------------------------------------------------------------------

figure;
subplot(1, 2, 1); hold on;
useCORAcolors("CORA:manual")
plot(S);

enlargeAxis(1.2)
title('$\mathcal{S}$','Interpreter','latex');
xlabel('$x_{(1)}$','Interpreter','latex')
ylabel('$x_{(2)}$','Interpreter','latex')

subplot(1, 2, 2); hold on;
useCORAcolors("CORA:manual-result")
plot(res)

enlargeAxis(1.2)
title('$M\otimes\mathcal{S}$','Interpreter','latex');
xlabel('$x_{(1)}$','Interpreter','latex')
% ylabel('$x_{(2)}$','Interpreter','latex')

% ------------------------------ END OF CODE ------------------------------
