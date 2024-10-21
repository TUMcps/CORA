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

% 1. comparison to empty set
% create polyZonotopes
c = [2; -3];
G = [2, 3; -1, 0];
GI = [1, 2, 4; 5, 6, 0];
E = [2, 0; 0, 1];

pZ1 = polyZonotope(c,G,GI,E);
assert(~representsa(pZ1,'emptySet'));

pZ2 = polyZonotope(c,G,[],E);
assert(~representsa(pZ2,'emptySet'));

pZ3 = polyZonotope.empty(2);
assert(representsa(pZ3,'emptySet'));


% 2. comparison to origin
% empty polyZonotope
pZ = polyZonotope.empty(2);
assert(~representsa(pZ,'origin'));

% only origin
pZ = polyZonotope(zeros(3,1));
assert(representsa(pZ,'origin'));

% shifted center
pZ = polyZonotope(0.01*ones(4,1));
assert(~representsa(pZ,'origin'));

% ...add tolerance
tol = 0.02;
assert(representsa(pZ,'origin',tol));

% include dependent generator matrix
pZ = polyZonotope(ones(2,1),0.1*eye(2));
tol = 2;
assert(representsa(pZ,'origin',tol));

% dependent and independent generators
pZ = polyZonotope(zeros(2,1),[0.01 -0.02; 0.03 0.01],[0.05; -0.02],[2 0; 0 1]);
tol = 0.1;
assert(representsa(pZ,'origin',tol));


% 3. comparison to zonotope
c = [2;-1]; G = [1;-1]; GI = [2 0 1; -1 1 0];

pZ = polyZonotope(c,[],GI);
[res,Z] = representsa(pZ,'zonotope');
assert(res)
assert(isequal(Z,zonotope(c,GI)));

pZ = polyZonotope(c,G,GI);
[res,Z] = representsa(pZ,'zonotope');
assert(res)
assert(isequal(Z,zonotope(c,[G,GI])));

E = [3; 0];
pZ = polyZonotope(c,G,GI,E);
[res,Z] = representsa(pZ,'zonotope');
assert(res)
assert(isequal(Z,zonotope(c,[G,GI])));


% 4. comparison to point
c = [2;-1];
G = [1;-1]; G_zeros = zeros(2,2);
GI = [2 0 1; -1 1 0]; GI_zeros = zeros(2,3);

pZ = polyZonotope(c);
assert(representsa(pZ,'point'));

pZ = polyZonotope(c,G);
assert(~representsa(pZ,'point'));

pZ = polyZonotope(c,G_zeros);
assert(representsa(pZ,'point'));

pZ = polyZonotope(c,[],GI);
assert(~representsa(pZ,'point'));

pZ = polyZonotope(c,[],GI_zeros);
assert(representsa(pZ,'point'));


% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
