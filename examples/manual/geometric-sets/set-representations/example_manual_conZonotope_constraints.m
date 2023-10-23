function example_manual_conZonotope_constraints()
% example_manual_conZonotope_constraints - example from the manual demonstrating the 
% constraints of a conZonotope as defined in the manual
%
% Syntax:
%   example_manual_conZonotope_constraints()
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
P = polytope(-A,-b);

% plot --------------------------------------------------------------------

figure; hold on; grid on;
xlim([-1;1]); ylim([-1;1]); zlim([-1;1]);
plot(P,1:3,'FaceColor','next');
view(31,21);

% enlargeAxis(1.2)
% title('$\mathcal{CZ}$','Interpreter','latex');
xlabel('$\beta_{1}$','Interpreter','latex')
ylabel('$\beta_{2}$','Interpreter','latex')
zlabel('$\beta_{3}$','Interpreter','latex')

% ------------------------------ END OF CODE ------------------------------
