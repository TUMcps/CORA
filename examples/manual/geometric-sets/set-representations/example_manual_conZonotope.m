function example_manual_conZonotope()
% example_manual_conZonotope - example from the manual demonstrating the 
% conZonotope constructor as defined in the manual
%
% Syntax:
%   example_manual_conZonotope()
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

% construct constrained zonotope
c = [0;0];
G = [1 0 1; 1 2 -1];
A = [-2 1 -1];
b = 2;

cZ = conZonotope(c,G,A,b);


% plot --------------------------------------------------------------------

figure; hold on
plot(cZ);

enlargeAxis(1.2)
title('$\mathcal{CZ}$','Interpreter','latex');
xlabel('$x_{(1)}$','Interpreter','latex')
ylabel('$x_{(2)}$','Interpreter','latex')

% ------------------------------ END OF CODE ------------------------------
