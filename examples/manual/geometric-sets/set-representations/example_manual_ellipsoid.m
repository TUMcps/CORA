function example_manual_ellipsoid()
% example_manual_ellipsoid - example from the manual demonstrating the 
% ellipsoid constructor as defined in the manual
%
% Syntax:
%   example_manual_ellipsoid()
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

% construct ellipsoid
Q = [13 7; 7 5];
q = [1; 2];

E = ellipsoid(Q,q);

% plot --------------------------------------------------------------------

figure; hold on
plot(E);

enlargeAxis(1.1)
title('$\mathcal{E}$','Interpreter','latex');
xlabel('$x_{(1)}$','Interpreter','latex')
ylabel('$x_{(2)}$','Interpreter','latex')

% ------------------------------ END OF CODE ------------------------------
