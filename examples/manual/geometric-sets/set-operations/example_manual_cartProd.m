function example_manual_cartProd()
% example_manual_cartProd - example from the manual demonstrating the
% cartProd operation as defined in the manual
%
% Syntax:
%   example_manual_cartProd()
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

% set S1 and S2
S1 = interval(-2,1);
S2 = interval(-1,2);
% Cartesian product
res = cartProd(S1,S2);

% plot --------------------------------------------------------------------

figure; hold on;
useCORAcolors("CORA:manual")
plot(res)

enlargeAxis(1.2)
title('$\mathcal{S}_1\times\mathcal{S}_2$','Interpreter','latex');
xlabel('$x_{(1)}$','Interpreter','latex')
ylabel('$x_{(2)}$','Interpreter','latex')

% ------------------------------ END OF CODE ------------------------------
