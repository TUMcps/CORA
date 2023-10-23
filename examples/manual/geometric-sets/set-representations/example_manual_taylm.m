function example_manual_taylm()
% example_manual_taylm - example from the manual demonstrating the 
% taylm constructor as defined in the manual
%
% Syntax:
%   example_manual_taylm()
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
f = @(x) cos(x);

% create Taylor model
D = interval(-1,1);
tx = taylm(D,2,'x');
ty = f(tx);

tay = [tx;ty];

% plot --------------------------------------------------------------------

figure; hold on
plot(tay,1:2,'FaceColor',CORAcolor('CORA:next'));

xs = linspace(-1,1,1000);
ys = f(xs);
plot(xs, ys, 'k');

xlabel('$x$','Interpreter','latex')
ylabel('$f(x)$','Interpreter','latex')

% ------------------------------ END OF CODE ------------------------------
