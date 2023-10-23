function example_manual_randPoint()
% example_manual_randPoint - example from the manual demonstrating the 
% randPoint operation as defined in the manual
%
% Syntax:
%   example_manual_randPoint()
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

% set seed
rng(1)

% set S
S = zonotope([0 1 1 0; 0 1 0 1]);

% random point
p = randPoint(S)

% plot --------------------------------------------------------------------

figure; hold on;
useCORAcolors("CORA:manual")
plot(S)
scatter(p(1,:),p(2,:),'.k')

enlargeAxis(1.2)
title('$\mathcal{S}$ and point','Interpreter','latex');
xlabel('$x_{(1)}$','Interpreter','latex')
ylabel('$x_{(2)}$','Interpreter','latex')

% ------------------------------ END OF CODE ------------------------------
