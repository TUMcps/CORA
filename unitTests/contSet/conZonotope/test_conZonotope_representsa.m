function res = test_conZonotope_representsa
% test_conZonotope_representsa - unit test function of representsa
%
% Syntax:
%    res = test_conZonotope_representsa
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
% Written:       21-April-2023
% Last update:   ---
% Last revision: 20-July-2023 (MW, rename '...representsa')

% ------------------------------ BEGIN CODE -------------------------------

% 1. comparison to empty set
% check empty conZonotope object
cZ = conZonotope.empty(2);
res = representsa(cZ,'emptySet');

% constrained zonotope
Z = [0 3 0 1;0 0 2 1];
A = [1 0 1]; b = 1;
cZ = conZonotope(Z,A,b);
res(end+1,1) = ~representsa(cZ,'emptySet');


% 2. comparsion to origin
% empty case
cZ = conZonotope.empty(2);
res(end+1,1) = ~representsa(cZ,'origin');

% true cases
cZ = conZonotope(zeros(3,1));
res(end+1) = representsa(cZ,'origin');

cZ = conZonotope(zeros(2,1),zeros(2));
res(end+1) = representsa(cZ,'origin');

% shifted center
cZ = conZonotope(ones(2,1),zeros(2,1));
res(end+1) = ~representsa(cZ,'origin');

% zero-centered, but with generator
cZ = conZonotope(zeros(3,1),[1; 0; -1]);
res(end+1) = ~representsa(cZ,'origin');

% below tolerance
cZ = conZonotope([0.5;0.5],[0.02 -0.1; 0.2 0.05]);
tol = 1;
res(end+1) = representsa_(cZ,'origin',tol);


% 3. comparison to zonotope
c = [1;2]; G = [1 -1 0; 2 1 3];
cZ = conZonotope(c,G);
[res(end+1),Z] = representsa(cZ,'zonotope');
res(end+1,1) = isequal(Z,zonotope(c,G));

Z = [0 3 0 1; 0 0 2 1];
A = [1 0 1]; b = 1;
cZ = conZonotope(Z,A,b);
res(end+1,1) = ~representsa(cZ,'zonotope');


% 4. comparison to interval
c = [1;2]; G = [1 0; 0 2];
cZ = conZonotope(c,G);
[res(end+1),I] = representsa(cZ,'interval');
res(end+1,1) = isequal(I,interval([0;0],[2;4]));


% 5. comparison to point
c = [1;2]; G = [1 0; 0 2];
cZ = conZonotope(c);
res(end+1,1) = representsa(cZ,'point');
cZ = conZonotope([c,G]);
res(end+1,1) = ~representsa(cZ,'point');


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
