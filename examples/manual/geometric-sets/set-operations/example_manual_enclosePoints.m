function example_manual_enclosePoints()
% example_manual_enclosePoints - example from the manual demonstrating the 
% enclosePoints operation as defined in the manual
%
% Syntax:
%   example_manual_enclosePoints()
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

% random point cloud
mu = [0 0];
sigma = [0.3 0.4; 0.4 1];
points = mvnrnd(mu,sigma,100)';

% compute enclosing set
S = ellipsoid.enclosePoints(points);

% plot --------------------------------------------------------------------

figure; hold on;
useCORAcolors("CORA:manual")
plot(S)
scatter(points(1,:),points(2,:),'.k')

enlargeAxis(1.2)
title('$\mathcal{S}$ and points','Interpreter','latex');
xlabel('$x_{(1)}$','Interpreter','latex')
ylabel('$x_{(2)}$','Interpreter','latex')

% ------------------------------ END OF CODE ------------------------------
