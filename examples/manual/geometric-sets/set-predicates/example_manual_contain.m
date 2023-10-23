function example_manual_contain()
% example_manual_contain - example from the manual demonstrating the 
% contains operation as defined in the manual
%
% Syntax:
%   example_manual_contain()
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

% sets S1,S2, and point p
S1 = zonotope([0 1 1 0; 0 1 0 1]);
S2 = interval([-1;-1],[1;1]);
p = [0.5;0.5];

% containment check
res1 = contains(S1,S2)
res2 = contains(S1,p)

% plot --------------------------------------------------------------------

figure; hold on;
useCORAcolors("CORA:manual")
plot(S1)
plot(S2)
scatter(p(1,:),p(2,:),'.k')

enlargeAxis(1.2)
title('$\mathcal{S}_1$, $\mathcal{S}_2$ and $p$','Interpreter','latex');
xlabel('$x_{(1)}$','Interpreter','latex')
ylabel('$x_{(2)}$','Interpreter','latex')

% ------------------------------ END OF CODE ------------------------------
