function res = test_contSet_full
% test_contSet_full - unit test function of full
%
% Syntax:
%    res = test_contSet_full
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
I = interval(sparse([1;2]),sparse([3;4]));
I_full = full(I);
assert(~issparse(I_full.inf));
assert(~issparse(I_full.sup));

% test zonotope
Z = zonotope(sparse([1;0]),sparse([1 0; 0 1]));
Z_full = full(Z);
assert(~issparse(Z_full.c));
assert(~issparse(Z_full.G));

% test capsule
C = capsule(sparse([1;0]),sparse([1;1]),sparse(0.5));
C_full = full(C);
assert(~issparse(C_full.c));
assert(~issparse(C_full.g));
assert(~issparse(C_full.r));

% test ellipsoid
E = ellipsoid(sparse([1 0; 0 2]),sparse([1;0]));
E_full = full(E);
assert(~issparse(E_full.Q));
assert(~issparse(E_full.q));

% test polytope (H-rep)
P = polytope(sparse([1 0; -1 0; 0 1; 0 -1]),sparse([1;1;1;1]));
P_full = full(P);
assert(~issparse(P_full.A));
assert(~issparse(P_full.b));

% test polytope (V-rep)
P = polytope(sparse([1 0; 0 1; -1 0]));
P_full = full(P);
assert(~issparse(P_full.V));

% test polyZonotope
pZ = polyZonotope(sparse([1;0]),sparse([1 0; 0 1]),sparse([0.5;0]),[1 0; 0 1]);
pZ_full = full(pZ);
assert(~issparse(pZ_full.c));
assert(~issparse(pZ_full.G));
assert(~issparse(pZ_full.GI));

% test conZonotope
cZ = conZonotope(sparse([0;0]),sparse([1 0; 0 1]),sparse([1 -1]),sparse(1));
cZ_full = full(cZ);
assert(~issparse(cZ_full.c));
assert(~issparse(cZ_full.G));

% test zonoBundle
Z1 = zonotope(sparse([1;0]),sparse([1 0; 0 1]));
Z2 = zonotope(sparse([0;0]),sparse([2 0; 0 2]));
zB = zonoBundle({Z1,Z2});
zB_full = full(zB);
assert(~issparse(zB_full.Z{1}.c));
assert(~issparse(zB_full.Z{2}.c));

% test emptySet
O = emptySet(2);
O_full = full(O);
assert(isa(O_full,'emptySet'));

% test fullspace
fs = fullspace(2);
fs_full = full(fs);
assert(isa(fs_full,'fullspace'));

% test that full of full keeps full
Z = zonotope([1;0],[1 0; 0 1]);
Z_full = full(Z);
assert(~issparse(Z_full.c));
assert(~issparse(Z_full.G));

% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
