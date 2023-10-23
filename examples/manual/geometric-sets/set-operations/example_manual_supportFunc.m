function example_manual_supportFunc()
% example_manual_supportFunc - example from the manual demonstrating the 
% supportFunc operation as defined in the manual
%
% Syntax:
%   example_manual_supportFunc()
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

% set S and vector l
S = zonotope([0 1 1 0; 0 1 0 1]);
l = [1;2];

% compute support function
res = supportFunc(S,l)

% plot --------------------------------------------------------------------

figure; hold on;
useCORAcolors("CORA:manual")
plot(S)
enlargeAxis(1.2)
plot(conHyperplane(-l,-res))

title('$\mathcal{S}$ and \texttt{supportFunc}','Interpreter','latex');
xlabel('$x_{(1)}$','Interpreter','latex')
ylabel('$x_{(2)}$','Interpreter','latex')

% ------------------------------ END OF CODE ------------------------------
