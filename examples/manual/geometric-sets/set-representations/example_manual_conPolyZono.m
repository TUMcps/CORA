function example_manual_conPolyZono()
% example_manual_conPolyZono - example from the manual demonstrating the 
% conPolyZono constructor as defined in the manual
%
% Syntax:
%   example_manual_conPolyZono()
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

% construct conPolyZono object
c = [0;0];
G = [1 0 1 -1;0 1 1 1];
E = [1 0 1 2;0 1 1 0;0 0 1 1];
A = [1 -0.5 0.5];
b = 0.5;
R = [0 1 2;1 0 0;0 1 0];

cPZ = conPolyZono(c,G,E,A,b,R);

% plot --------------------------------------------------------------------

figure; hold on
plot(cPZ);

enlargeAxis(1.1)
title('$\mathcal{CPZ}$','Interpreter','latex');
xlabel('$x_{(1)}$','Interpreter','latex')
ylabel('$x_{(2)}$','Interpreter','latex')

% ------------------------------ END OF CODE ------------------------------
