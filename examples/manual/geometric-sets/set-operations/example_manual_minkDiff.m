function example_manual_minkDiff()
% example_manual_minkDiff - example from the manual demonstrating the
% minkDiff operation as defined in the manual
%
% Syntax:
%   example_manual_minkDiff()
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

% set 1 and set 2
S1 = zonotope([0.5 0.5 -0.3 1 0; 0 0.2 1 0 1]);
S2 = interval([-1;-1],[1;1]);

% Minkowski difference
res = minkDiff(S1,S2);

% plot --------------------------------------------------------------------

figure;
subplot(1, 2, 1); hold on;
useCORAcolors("CORA:manual")
plot(S1);
plot(S2);

title('$\mathcal{S}_1$ and $\mathcal{S}_2$','Interpreter','latex');
xlabel('$x_{(1)}$','Interpreter','latex')
ylabel('$x_{(2)}$','Interpreter','latex')
xlim([-3,3]);ylim([-3,3]);

subplot(1, 2, 2); hold on;
useCORAcolors("CORA:manual-result")
plot(res)

title('$\texttt{minkDiff}(\mathcal{S}_1,\mathcal{S}_2)$','Interpreter','latex');
xlabel('$x_{(1)}$','Interpreter','latex')
% ylabel('$x_{(2)}$','Interpreter','latex')
xlim([-3,3]);ylim([-3,3]);

% ------------------------------ END OF CODE ------------------------------
