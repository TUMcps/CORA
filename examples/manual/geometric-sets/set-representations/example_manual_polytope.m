function example_manual_polytope()
% example_manual_polytope - example from the manual demonstrating the 
% polytope constructor as defined in the manual
%
% Syntax:
%   example_manual_polytope()
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

 % construct polytope (halfspace rep.)
 C = [1 0 -1 0 1; 0 1 0 -1 1]';
 d = [3; 2; 3; 2; 1];
 poly = polytope(C,d);

 % construct polytope (vertex rep.)
 V = [-3 -3 -1 3; -2 2 2 -2];
 poly = polytope(V);

% plot --------------------------------------------------------------------

figure; hold on
plot(poly);

enlargeAxis(1.2)
title('$\mathcal{P}$','Interpreter','latex');
xlabel('$x_{(1)}$','Interpreter','latex')
ylabel('$x_{(2)}$','Interpreter','latex')

% ------------------------------ END OF CODE ------------------------------
