function example_manual_example_conZonotope_zonotope()
% example_manual_example_conZonotope_zonotope - example from the manual
% demonstrating the zonotope of a conZonotope as defined in the manual
%
% Syntax:
%   example_manual_example_conZonotope_zonotope()
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

% construct constrained zonotope
c = [0;0];
G = [1 0 1; 1 2 -1];
A = [-2 1 -1];
b = 2;

cZ = conZonotope(c,G,A,b);
Z = zonotope(c,G);

% plot --------------------------------------------------------------------

figure; hold on
plot(cZ);
plot(Z);

enlargeAxis(1.2)
% title('$\mathcal{CZ}$','Interpreter','latex');
xlabel('$x_{(1)}$','Interpreter','latex')
ylabel('$x_{(2)}$','Interpreter','latex')

% ------------------------------ END OF CODE ------------------------------
