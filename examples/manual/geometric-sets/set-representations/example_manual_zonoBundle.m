function example_manual_zonoBundle()
% example_manual_zonoBundle - example from the manual demonstrating the 
% zonoBundle constructor as defined in the manual
%
% Syntax:
%   example_manual_zonoBundle()
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

% construct zonotopes
zono1 = zonotope([1 3 0; 1 0 2]);
zono2 = zonotope([0 2 2; 0 2 -2]);

% construct zonotope bundle
list = {zono1,zono2};

zB = zonoBundle(list);

% plot --------------------------------------------------------------------

figure; hold on
plot(zB);

enlargeAxis(1.2)
title('$\mathcal{ZB}$','Interpreter','latex');
xlabel('$x_{(1)}$','Interpreter','latex')
ylabel('$x_{(2)}$','Interpreter','latex')

% ------------------------------ END OF CODE ------------------------------
