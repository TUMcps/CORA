function example_manual_center()
% example_manual_center - example from the manual demonstrating the 
% center operation as defined in the manual
%
% Syntax:
%   example_manual_center()
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

% set S
S = interval([-2;-2],[1;1]);

% compute center
res = center(S)

% plot --------------------------------------------------------------------

figure; hold on;
useCORAcolors("CORA:manual")
plot(S)
scatter(res(1,:),res(2,:),'.k')

enlargeAxis(1.2)
title('$\mathcal{S}$ and center','Interpreter','latex');
xlabel('$x_{(1)}$','Interpreter','latex')
ylabel('$x_{(2)}$','Interpreter','latex')

% ------------------------------ END OF CODE ------------------------------
