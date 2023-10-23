function example_manual_zonotope()
% example_manual_zonotope - example from the manual demonstrating the 
% zonotope constructor as defined in the manual
%
% Syntax:
%   example_manual_zonotope()
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

% construct zonotope
c = [1;1];
G = [1 1 1; 1 -1 0];

zono = zonotope(c,G);

% plot --------------------------------------------------------------------

figure; hold on
plot(zono);

enlargeAxis(1.2)
title('$\mathcal{Z}$','Interpreter','latex');
xlabel('$x_{(1)}$','Interpreter','latex')
ylabel('$x_{(2)}$','Interpreter','latex')

% ------------------------------ END OF CODE ------------------------------
