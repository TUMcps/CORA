function example_manual_plus()
% example_manual_plus - example from the manual demonstrating the
% plus operation as defined in the manual
%
% Syntax:
%   example_manual_plus()
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

% set S1 and S2
S1 = zonotope([0 0.5 1; 0 1 0]);
S2 = zonotope([0 1 0; 0 0 1]);

% Minkowski sum
res = S1 + S2;

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

title('$\mathcal{S}_1\oplus\mathcal{S}_2$','Interpreter','latex');
xlabel('$x_{(1)}$','Interpreter','latex')
% ylabel('$x_{(2)}$','Interpreter','latex')
xlim([-3,3]);ylim([-3,3]);

% ------------------------------ END OF CODE ------------------------------
