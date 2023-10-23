function example_manual_example_conPolyZono()
% example_manual_example_conPolyZono - example from the manual demonstrating 
% the conPolyZono example from the manual
%
% Syntax:
%   example_manual_example_conPolyZono()
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
f1 = figure; hold on
plot(cPZ,[1,2],'b');

f2 = figure; hold on
plot(res,[1,2],'r','Splits',25);

% plot --------------------------------------------------------------------

figure(f1)
enlargeAxis(1.1)
xlabel('$x_{(1)}$','Interpreter','latex')
ylabel('$x_{(2)}$','Interpreter','latex')

figure(f2)
enlargeAxis(1.1)
xlabel('$x_{(1)}$','Interpreter','latex')
ylabel('$x_{(2)}$','Interpreter','latex')

% ------------------------------ END OF CODE ------------------------------
