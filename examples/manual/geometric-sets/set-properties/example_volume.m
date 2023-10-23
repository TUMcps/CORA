% set S
S = zonotope([0 1 1 0; 0 1 0 1]);

% volume of the set
res = volume(S)

% plot --------------------------------------------------------------------

figure; hold on;
useCORAcolors("CORA:manual")
plot(S)

enlargeAxis(1.2)
title('$\mathcal{S}$','Interpreter','latex');
xlabel('$x_{(1)}$','Interpreter','latex')
ylabel('$x_{(2)}$','Interpreter','latex')