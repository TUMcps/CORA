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

% empty conZonotope
n = 2;
cZ = conZonotope.empty(n);
assert(representsa_(cZ,'emptySet',eps) && dim(cZ) == n);
cZ = conZonotope(zeros(3,0));
assert(representsa_(cZ,'emptySet',eps))
assert(size(cZ.c,1) == 3 && size(cZ.G,1) == 3);


% init simple constrained zonotope
Z = [0 3 0 1;0 0 2 1];
A = [1 0 1];
b = 1;
cZ = conZonotope(Z,A,b);

assert(all(withinTol([cZ.c,cZ.G],Z),'all') ...
        && all(withinTol(cZ.A,A),'all') ...
        && all(withinTol(cZ.b,b)));

% constraints should be double for internal processing
c = [ 1.000 ; 1.000 ];
G = [ 3.000 0.000 0.000 ; 0.000 3.000 0.000 ];
A = single([ 1.211 -0.118 0.684 ]);
b = single(-0.112);
cZ = conZonotope(c,G,A,b);
assert(isa(cZ.A,'double') && isa(cZ.b,'double'))


% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
