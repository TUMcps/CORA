function example_manual_probZonotope()
% example_manual_probZonotope - example from the manual demonstrating the 
% probZonotope constructor as defined in the manual
%
% Syntax:
%   example_manual_probZonotope()
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

% construct probabilistic zonotope
c = [0;0];
G = [1 0;0 1];
G_ = [3 2; 3 -2];

probZ = probZonotope([c,G],G_);

% plot --------------------------------------------------------------------

figure; hold on; grid on;
plot(probZ);
view(31,21);

enlargeAxis(1.2)
title('$\mathscr{Z}$','Interpreter','latex');
xlabel('$x_{(1)}$','Interpreter','latex')
ylabel('$x_{(2)}$','Interpreter','latex')

% ------------------------------ END OF CODE ------------------------------
