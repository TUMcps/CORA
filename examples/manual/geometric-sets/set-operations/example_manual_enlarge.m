function example_manual_enlarge()
% example_manual_enlarge - example from the manual demonstrating the 
% enlarge operation as defined in the manual
%
% Syntax:
%   example_manual_enlarge()
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
% Written:       05-February-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% set S and scaling lambda
S = zonotope([1;1],[1 0 1; 0 1 1]);
lambda = 1.2;

% enlarge set S
S_ = enlarge(S,lambda);

% plot --------------------------------------------------------------------

figure; hold on;
useCORAcolors("CORA:manual")
plot(S)
plot(S_)

enlargeAxis(1.2)
title('$\mathcal{S}$ and $\mathcal{S}''$','Interpreter','latex');
xlabel('$x_{(1)}$','Interpreter','latex')
ylabel('$x_{(2)}$','Interpreter','latex')

% ------------------------------ END OF CODE ------------------------------
