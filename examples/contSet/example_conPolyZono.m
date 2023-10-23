function completed = example_conPolyZono()
% example_conPolyZono - example demonstrating set based computation with 
%                       constrained polynomial zonotopes
%
% Syntax:
%    completed = example_conPolyZono()
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
% See also: example_polyZonotope

% Authors:       Niklas Kochdumper
% Written:       08-February-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% construct zonotope
Z = zonotope([0;0],[1 1;0 1]);

% construct ellipsoid
E = ellipsoid([2 1;1 2],[1;1]);

% convert sets to constrained polynomial zonotopes
cPZ1 = conPolyZono(Z);
cPZ2 = conPolyZono(E);

% compute the Minkowski sum
resSum = cPZ1 + cPZ2;

% compute the intersection
resAnd = cPZ1 & cPZ2;

% compute the union
resOR = cPZ1 | cPZ2;

% construct conPolyZono object
c = [0;0];
G = [1 0 1 -1;0 1 1 1];
E = [1 0 1 2;0 1 1 0;0 0 1 1];
A = [1 -0.5 0.5];
b = 0.5;
R = [0 1 2;1 0 0;0 1 0];
 
cPZ = conPolyZono(c,G,E,A,b,R);

% compute quadratic map
Q{1} = [0.5 0.5; 0 -0.5];
Q{2} = [-1 0; 1 1];

res = quadMap(cPZ,Q);

% visualization
figure;
subplot(1,2,1); hold on;
plot(cPZ);

subplot(1,2,2); hold on;
plot(res,[1,2],'Color',colorblind('r'),'Splits',25);

% example completed
completed = true;

% ------------------------------ END OF CODE ------------------------------
