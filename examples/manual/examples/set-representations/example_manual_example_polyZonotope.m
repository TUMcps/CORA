function example_manual_example_polyZonotope()
% example_manual_example_polyZonotope - example from the manual demonstrating 
% the polyZonotope example from the manual
%
% Syntax:
%   example_manual_example_polyZonotope()
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

% construct zonotope
c = [1;0];
G = [1 1;1 0];
Z = zonotope(c,G);

% compute over-approximation of the quadratic map
Q{1} = [0.5 0.5; 0 -0.5];
Q{2} = [-1 0; 1 1];
resZono = quadMap(Z,Q);

% convert zonotope to polynomial zonotope
pZ = polyZonotope(Z);

% compute the exact quadratic map
resPolyZono = quadMap(pZ,Q);

% visualization
figure; hold on;
plot(resZono,[1,2],'r');
plot(resPolyZono,[1,2],'b');

% plot --------------------------------------------------------------------

enlargeAxis(1.2)
xlabel('$x_{(1)}$','Interpreter','latex')
ylabel('$x_{(2)}$','Interpreter','latex')

% ------------------------------ END OF CODE ------------------------------
