function example_manual_zonotope_construction()
% example_manual_zonotope_construction - example from the manual
% demonstrating the construction of a zonotope as defined in the manual
%
% Syntax:
%   example_manual_zonotope_construction()
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

Z1 = zonotope(c,G(:,1:1));
Z2 = zonotope(c,G(:,1:2));
Z3 = zonotope(c,G(:,1:3));

% plot --------------------------------------------------------------------

figure; 

subplot(1, 3, 1); hold on
plot(Z1);
V = [c vertices(Z1)];
scatter(V(1,:),V(2,:),'.k')
text(1.2,1,'$c$','Interpreter','latex')
text(1.4,1.6,'$l^{(1)}$','Interpreter','latex')

enlargeAxis(1.2)
title('$c\oplus l^{(1)}$','Interpreter','latex');
xlabel('$x_{(1)}$','Interpreter','latex')
ylabel('$x_{(2)}$','Interpreter','latex')

subplot(1, 3, 2); hold on
plot(Z2);
V = [c vertices(Z2)];
scatter(V(1,:),V(2,:),'.k')
text(1.3,1,'$c$','Interpreter','latex')
text(-0.5,2,'$l^{(1)}$','Interpreter','latex')
text(1.4,2,'$l^{(2)}$','Interpreter','latex')

enlargeAxis(1.2)
title('$c\oplus l^{(1)}\oplus l^{(2)}$','Interpreter','latex');
xlabel('$x_{(1)}$','Interpreter','latex')
% ylabel('$x_{(2)}$','Interpreter','latex')

subplot(1, 3, 3); hold on
plot(Z3);
V = [c vertices(Z3)];
scatter(V(1,:),V(2,:),'.k')
text(1.4,1,'$c$','Interpreter','latex')
text(-1.4,2,'$l^{(1)}$','Interpreter','latex')
text(2.2,2,'$l^{(2)}$','Interpreter','latex')
text(1,2.7,'$l^{(3)}$','Interpreter','latex')

enlargeAxis(1.2)
title('$c\oplus l^{(1)}\oplus l^{(2)}\oplus l^{(3)}$','Interpreter','latex');
xlabel('$x_{(1)}$','Interpreter','latex')
% ylabel('$x_{(2)}$','Interpreter','latex')

% ------------------------------ END OF CODE ------------------------------
