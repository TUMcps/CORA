function res = test_contSet_double
% test_contSet_double - unit test function of double
%
% Syntax:
%    res = test_contSet_double
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
% See also: none

% Authors:       Tobias Ladner
% Written:       02-April-2026
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% test interval
I = interval(single([1;2]),single([3;4]));
I_dbl = double(I);
assert(isa(I_dbl.inf,'double'));
assert(isa(I_dbl.sup,'double'));

% test zonotope
Z = zonotope(single([1;0]),single([1 0; 0 1]));
Z_dbl = double(Z);
assert(isa(Z_dbl.c,'double'));
assert(isa(Z_dbl.G,'double'));

% test capsule
C = capsule(single([1;0]),single([1;1]),single(0.5));
C_dbl = double(C);
assert(isa(C_dbl.c,'double'));
assert(isa(C_dbl.g,'double'));
assert(isa(C_dbl.r,'double'));

% test ellipsoid
E = ellipsoid(single([1 0; 0 2]),single([1;0]));
E_dbl = double(E);
assert(isa(E_dbl.Q,'double'));
assert(isa(E_dbl.q,'double'));

% test polytope (H-rep)
P = polytope(single([1 0; -1 0; 0 1; 0 -1]),single([1;1;1;1]));
P_dbl = double(P);
assert(isa(P_dbl.A,'double'));
assert(isa(P_dbl.b,'double'));

% test polytope (V-rep)
P = polytope(single([1 0; 0 1; -1 0]));
P_dbl = double(P);
assert(isa(P_dbl.V,'double'));

% test polyZonotope
pZ = polyZonotope(single([1;0]),single([1 0; 0 1]),single([0.5;0]),[1 0; 0 1]);
pZ_dbl = double(pZ);
assert(isa(pZ_dbl.c,'double'));
assert(isa(pZ_dbl.G,'double'));
assert(isa(pZ_dbl.GI,'double'));

% test conZonotope
cZ = conZonotope(single([0;0]),single([1 0; 0 1]),single([1 -1]),single(1));
cZ_dbl = double(cZ);
assert(isa(cZ_dbl.c,'double'));
assert(isa(cZ_dbl.G,'double'));

% test zonoBundle
Z1 = zonotope(single([1;0]),single([1 0; 0 1]));
Z2 = zonotope(single([0;0]),single([2 0; 0 2]));
zB = zonoBundle({Z1,Z2});
zB_dbl = double(zB);
assert(isa(zB_dbl.Z{1}.c,'double'));
assert(isa(zB_dbl.Z{2}.c,'double'));

% test emptySet
O = emptySet(2);
O_dbl = double(O);
assert(isa(O_dbl,'emptySet'));

% test fullspace
fs = fullspace(2);
fs_dbl = double(fs);
assert(isa(fs_dbl,'fullspace'));

% test that double of double keeps double
Z = zonotope([1;0],[1 0; 0 1]);
Z_dbl = double(Z);
assert(isa(Z_dbl.c,'double'));
assert(isa(Z_dbl.G,'double'));

% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
