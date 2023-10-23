function example_manual_polyZonotope()
% example_manual_polyZonotope - example from the manual demonstrating the 
% polyZonotope constructor as defined in the manual
%
% Syntax:
%   example_manual_polyZonotope()
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

% construct polynomial zonotope
c = [4;4];
G = [2 1 2; 0 2 2];
E = [1 0 3;0 1 1];
GI = [1;0];

pZ = polyZonotope(c,G,GI,E);

% plot --------------------------------------------------------------------

figure; hold on
plot(pZ);

enlargeAxis(1.2)
title('$\mathcal{PZ}$','Interpreter','latex');
xlabel('$x_{(1)}$','Interpreter','latex')
ylabel('$x_{(2)}$','Interpreter','latex')

% ------------------------------ END OF CODE ------------------------------
