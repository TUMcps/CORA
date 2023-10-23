function example_manual_isIntersecting()
% example_manual_isIntersecting - example from the manual demonstrating the 
% isIntersecting operation as defined in the manual
%
% Syntax:
%   example_manual_isIntersecting()
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

% sets S1 and S2
S1 = interval([-1;-1],[2;2]);
S2 = interval([-2;-2],[1;1]);

% intersection check
res = isIntersecting(S1,S2)

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
