function example_manual_contract()
% example_manual_contract - example from the manual demontrating the
% contract operation as defined in the manual
%
% Syntax:
%   example_manual_contract()
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

% function f(x)
f = @(x) x(1)^2 + x(2)^2 - 4;

% domain D
dom = interval([1;1],[2.5;2.5]);

% contraction
res = contract(f,dom,'forwardBackward');

% plot --------------------------------------------------------------------

syms x y;
eqs = f([x;y]);
ls = levelSet(eqs,[x;y],'==');

figure; hold on;
xlim([-3,3]); ylim([-3,3]);

useCORAcolors("CORA:manual")
plot(ls)
plot(dom)

useCORAcolors("CORA:manual-result");
updateColorIndex(0);
plot(res)

% title('$\mathcal{S}$ and center','Interpreter','latex');
xlabel('$x_{(1)}$','Interpreter','latex')
ylabel('$x_{(2)}$','Interpreter','latex')

% ------------------------------ END OF CODE ------------------------------
