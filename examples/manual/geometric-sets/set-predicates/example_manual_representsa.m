function example_manual_representsa()
% example_manual_representsa - example from the manual demonstrating the 
% representsa operation as defined in the manual
%
% Syntax:
%   example_manual_representsa()
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

% set S (intersection)
S1 = polytope([-1 -1;0 -1;0 1;1 1], [-0.5; 0; 2; 2.5]);
S2 = polytope([-1 -1;0 -1;0 1;1 1], [2.5; 2; 0; -0.5]);
S = S1 & S2;

% check if set is empty
res = representsa(S,'emptySet')

% plot --------------------------------------------------------------------

figure; hold on;
useCORAcolors("CORA:manual")
plot(S1)
plot(S2)

enlargeAxis(1.2)
title('$\mathcal{S}_1$ and $\mathcal{S}_2$','Interpreter','latex');
xlabel('$x_{(1)}$','Interpreter','latex')
ylabel('$x_{(2)}$','Interpreter','latex')

% ------------------------------ END OF CODE ------------------------------
