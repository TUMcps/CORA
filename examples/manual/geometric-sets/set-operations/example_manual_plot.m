function example_manual_plot()
% example_manual_plot - example from the manual demonstrating the
% plot operation as defined in the manual
%
% Syntax:
%   example_manual_plot()
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
S = zonotope([0 1 1 0; 0 1 2 1; 0 1 0 1]);

% visualization
figure;
plot(S,[1,3],'--r');

% plot --------------------------------------------------------------------

enlargeAxis(1.2)
title('$\mathcal{S}$','Interpreter','latex');
xlabel('$x_{(1)}$','Interpreter','latex')
ylabel('$x_{(2)}$','Interpreter','latex')

% ------------------------------ END OF CODE ------------------------------
