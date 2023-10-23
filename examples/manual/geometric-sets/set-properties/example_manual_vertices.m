function example_manual_vertices()
% example_manual_vertices - example from the manual demonstrating the 
% vertices operation as defined in the manual
%
% Syntax:
%   example_manual_vertices()
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

% set S
S = interval([-2;-2], [1;1]);

% compute vertices
V = vertices(S)

% plot --------------------------------------------------------------------

figure; hold on;
useCORAcolors("CORA:manual")
plot(S)
scatter(V(1,:),V(2,:),'.k')

title('$\mathcal{S}$ and $V$','Interpreter','latex');
xlabel('$x_{(1)}$','Interpreter','latex')
ylabel('$x_{(2)}$','Interpreter','latex')
xlim([-3,2]); ylim([-3,2]);

% ------------------------------ END OF CODE ------------------------------
