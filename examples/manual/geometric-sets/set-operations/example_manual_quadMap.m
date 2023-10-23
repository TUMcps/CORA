function example_manual_quadMap()
% example_manual_quadMap - example from the manual demonstrating the
% quadMap operation as defined in the manual
%
% Syntax:
%   example_manual_quadMap()
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
Q{1} = [0.5 0.5; 0 -0.5];
Q{2} = [-1 0; 1 1];

% quadratic map
res = quadMap(S,Q);

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

title('$\texttt{quadMap}(\mathcal{S},Q)$','Interpreter','latex');
xlabel('$x_{(1)}$','Interpreter','latex')
% ylabel('$x_{(2)}$','Interpreter','latex')
xlim([-3,3]);ylim([-3,3]);

% ------------------------------ END OF CODE ------------------------------
