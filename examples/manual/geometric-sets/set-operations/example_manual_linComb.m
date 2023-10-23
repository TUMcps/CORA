function example_manual_linComb()
% example_manual_linComb - example from the manual demonstrating the 
% linComb operation as defined in the manual
%
% Syntax:
%   example_manual_linComb()
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
S1 = polyZonotope([0.5;0.5],[1 1;-1 1],[],[1 2]);
S2 = zonotope([-1.5;-1.5],[1 0;0 1]);

% linear combination
res = linComb(S1,S2);

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
plot(res,1:2,'Splits',20)

title('$\texttt{linComb}(\mathcal{S}_1,\mathcal{S}_2)$','Interpreter','latex');
xlabel('$x_{(1)}$','Interpreter','latex')
% ylabel('$x_{(2)}$','Interpreter','latex')
xlim([-3,3]);ylim([-3,3]);

% ------------------------------ END OF CODE ------------------------------
