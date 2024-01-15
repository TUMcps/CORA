function res = test_polyZonotope_representsa
% test_polyZonotope_representsa - unit test function of representsa
%
% Syntax:
%    res = test_polyZonotope_representsa
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Mark Wetzlinger
% Written:       01-May-2020
% Last update:   ---
% Last revision: 20-July-2023 (MW, rename '...representsa')

% ------------------------------ BEGIN CODE -------------------------------

res = true(0);

% 1. comparison to empty set
% create polyZonotopes
c = [2; -3];
G = [2, 3; -1, 0];
GI = [1, 2, 4; 5, 6, 0];
E = [2, 0; 0, 1];

pZ1 = polyZonotope(c,G,GI,E);
res(end+1,1) = ~representsa(pZ1,'emptySet');

pZ2 = polyZonotope(c,G,[],E);
res(end+1,1) = ~representsa(pZ2,'emptySet');

pZ3 = polyZonotope.empty(2);
res(end+1,1) = representsa(pZ3,'emptySet');


% 2. comparison to origin
% empty polyZonotope
pZ = polyZonotope.empty(2);
res(end+1,1) = ~representsa(pZ,'origin');

% only origin
pZ = polyZonotope(zeros(3,1));
res(end+1,1) = representsa(pZ,'origin');

% shifted center
pZ = polyZonotope(0.01*ones(4,1));
res(end+1,1) = ~representsa(pZ,'origin');

% ...add tolerance
tol = 0.02;
res(end+1,1) = representsa(pZ,'origin',tol);

% include dependent generator matrix
pZ = polyZonotope(ones(2,1),0.1*eye(2));
tol = 2;
res(end+1,1) = representsa(pZ,'origin',tol);

% dependent and independent generators
pZ = polyZonotope(zeros(2,1),[0.01 -0.02; 0.03 0.01],[0.05; -0.02],[2 0; 0 1]);
tol = 0.1;
res(end+1,1) = representsa(pZ,'origin',tol);


% 3. comparison to zonotope
c = [2;-1]; G = [1;-1]; GI = [2 0 1; -1 1 0];

pZ = polyZonotope(c,[],GI);
[res(end+1,1),Z] = representsa(pZ,'zonotope');
res(end+1,1) = isequal(Z,zonotope(c,GI));

pZ = polyZonotope(c,G,GI);
[res(end+1,1),Z] = representsa(pZ,'zonotope');
res(end+1,1) = isequal(Z,zonotope(c,[G,GI]));

E = [3; 0];
pZ = polyZonotope(c,G,GI,E);
[res(end+1,1),Z] = representsa(pZ,'zonotope');
res(end+1,1) = isequal(Z,zonotope(c,[G,GI]));


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
