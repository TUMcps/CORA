function example_manual_matrixSet_mtimes()
% example_manual_matrixSet_mtimes - example from the manual demonstrating 
% the mtimes operation of a matrix set as defined in the manual
%
% Syntax:
%   example_manual_matrixSet_mtimes()
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

% vector set
S = zonotope([0 1 1 0; 0 1 0 1]);

% matrix set
C = [1 0; -1 0.5];
D = [0.1 0; 0 0.2];
A = intervalMatrix(C,D);

% linear transformation
res = A*S;

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

title('$A\otimes\mathcal{S}$','Interpreter','latex');
xlabel('$x_{(1)}$','Interpreter','latex')
% ylabel('$x_{(2)}$','Interpreter','latex')
xlim([-3,3]);ylim([-3,3]);

% ------------------------------ END OF CODE ------------------------------
