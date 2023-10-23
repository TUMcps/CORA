function example_manual_capsule()
% example_manual_capsule - example from the manual demonstrating the 
% capsule constructor as defined in the manual
%
% Syntax:
%   example_manual_capsule()
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

% construct capsule
c = [1;2];
g = [2;1];
r = 1;

C = capsule(c,g,r);

% plot --------------------------------------------------------------------

figure; hold on
plot(C);

enlargeAxis(1.1)
title('$\mathcal{C}$','Interpreter','latex');
xlabel('$x_{(1)}$','Interpreter','latex')
ylabel('$x_{(2)}$','Interpreter','latex')

% ------------------------------ END OF CODE ------------------------------
