function res = test_conZonotope_conZonotope
% test_conZonotope_conZonotope - unit test function of conZonotope (constructor)
%
% Syntax:
%    res = test_conZonotope_conZonotope
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
% Written:       19-March-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true(0);

% empty conZonotope
n = 2;
cZ = conZonotope.empty(n);
res(end+1,1) = representsa_(cZ,'emptySet',eps) && dim(cZ) == n;
cZ = conZonotope(zeros(3,0));
res(end+1,1) = representsa_(cZ,'emptySet',eps) ...
    && size(cZ.c,1) == 3 && size(cZ.G,1) == 3;


% init simple constrained zonotope
Z = [0 3 0 1;0 0 2 1];
A = [1 0 1];
b = 1;
cZ = conZonotope(Z,A,b);

res(end+1,1) = all(all(withinTol([cZ.c,cZ.G],Z))) ...
        && all(all(withinTol(cZ.A,A))) ...
        && all(withinTol(cZ.b,b));


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
