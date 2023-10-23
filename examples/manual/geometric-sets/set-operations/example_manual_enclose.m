function example_manual_enclose()
% example_manual_enclose - example from the manual demonstrating the 
% enclose operation as defined in the manual
%
% Syntax:
%   example_manual_enclose()
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

% sets S1,S2 and matrix M
S1 = polyZonotope([1.5;1.5], ...
    [1 0;0 1], ...
    [],eye(2));
S2 = [0.5;0.5];
M = [-1 0;0 -1];

% apply method enclose
S3 = M*S1 + S2;

res = enclose(S1,M,S2);
res = enclose(S1,S3);

% plot --------------------------------------------------------------------

figure;
subplot(1, 2, 1); hold on;
useCORAcolors("CORA:manual")
plot(S1);
plot(S3);

title('$\mathcal{S}_1$ and $\mathcal{S}_3$','Interpreter','latex');
xlabel('$x_{(1)}$','Interpreter','latex')
ylabel('$x_{(2)}$','Interpreter','latex')
xlim([-3,3]);ylim([-3,3]);

subplot(1, 2, 2); hold on;
useCORAcolors("CORA:manual-result")
plot(res,1:2,'Splits',20)

title('$\texttt{enclose}(\mathcal{S}_1, \mathcal{S}_3)$','Interpreter','latex');
xlabel('$x_{(1)}$','Interpreter','latex')
% ylabel('$x_{(2)}$','Interpreter','latex')
xlim([-3,3]);ylim([-3,3]);

% ------------------------------ END OF CODE ------------------------------
