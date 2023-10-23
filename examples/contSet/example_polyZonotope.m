function completed = example_polyZonotope()
% example_polyZonotope - example demonstrating set based computation with 
%                        polynomial zonotopes
%
% Syntax:
%    completed = example_polyZonotope()
%
% Inputs:
%    -
%
% Outputs:
%    completed - true/false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       ---
% Written:       ---
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% construct zonotope
c = [1;0];
G = [1 1;1 0];
zono = zonotope(c,G);

% compute over-approximation of the quadratic map
Q{1} = [0.5 0.5; 0 -0.5];
Q{2} = [-1 0; 1 1];

resZono = quadMap(zono,Q);

% convert zonotope to polynomial zonotope
pZ = polyZonotope(zono);

% compute the exact quadratic map
resPolyZono = quadMap(pZ,Q);

% visualization
figure; hold on;
plot(resZono,[1,2],'Color',colorblind('r'));
plot(resPolyZono);

% example completed
completed = true;

% ------------------------------ END OF CODE ------------------------------
