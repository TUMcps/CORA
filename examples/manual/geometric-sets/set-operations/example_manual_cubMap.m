function example_manual_cubMap()
% example_manual_cubMap - example from the manual demonstrating the 
% cubMap operation as defined in the manual
%
% Syntax:
%   example_manual_cubMap()
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

% set and matrices
S = polyZonotope([0;0], [1 1;1 0], [],eye(2));

T{1,1} = 0.4*[1 2; -1 2];
T{1,2} = 0.4*[-3 0; 1 1];
T{2,1} = 0.05*[2 0; -2 1];
T{2,2} = 0.05*[-3 0; -21 -1];

% cubic map
res = cubMap(S,T);

% plot --------------------------------------------------------------------

figure;
subplot(1, 2, 1); hold on;
useCORAcolors("CORA:manual")
plot(S);

title('$\mathcal{S}$','Interpreter','latex');
xlabel('$x_{(1)}$','Interpreter','latex')
ylabel('$x_{(2)}$','Interpreter','latex')
xlim([-3,3]);ylim([-3,3]);

subplot(1, 2, 2); hold on;
useCORAcolors("CORA:manual-result")
plot(res)

title('$\texttt{cubMap}(\mathcal{S},T)$','Interpreter','latex');
xlabel('$x_{(1)}$','Interpreter','latex')
% ylabel('$x_{(2)}$','Interpreter','latex')
xlim([-3,3]);ylim([-3,3]);

% ------------------------------ END OF CODE ------------------------------
