function res = test_contSet_single
% test_contSet_single - unit test function of single
%
% Syntax:
%    res = test_contSet_single
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
I = interval([1;2],[3;4]);
I_sgl = single(I);
assert(isa(I_sgl.inf,'single'));
assert(isa(I_sgl.sup,'single'));

% test zonotope
Z = zonotope([1;0],[1 0; 0 1]);
Z_sgl = single(Z);
assert(isa(Z_sgl.c,'single'));
assert(isa(Z_sgl.G,'single'));

% test capsule
C = capsule([1;0],[1;1],0.5);
C_sgl = single(C);
assert(isa(C_sgl.c,'single'));
assert(isa(C_sgl.g,'single'));
assert(isa(C_sgl.r,'single'));

% test ellipsoid
E = ellipsoid([1 0; 0 2],[1;0]);
E_sgl = single(E);
assert(isa(E_sgl.Q,'single'));
assert(isa(E_sgl.q,'single'));

% test polytope (H-rep)
P = polytope([1 0; -1 0; 0 1; 0 -1],[1;1;1;1]);
P_sgl = single(P);
assert(isa(P_sgl.A,'single'));
assert(isa(P_sgl.b,'single'));

% test polytope (V-rep)
P = polytope([1 0; 0 1; -1 0]);
P_sgl = single(P);
assert(isa(P_sgl.V,'single'));

% test polyZonotope
pZ = polyZonotope([1;0],[1 0; 0 1],[0.5;0],[1 0; 0 1]);
pZ_sgl = single(pZ);
assert(isa(pZ_sgl.c,'single'));
assert(isa(pZ_sgl.G,'single'));
assert(isa(pZ_sgl.GI,'single'));

% test conZonotope
cZ = conZonotope([0;0],[1 0; 0 1],[1 -1],1);
cZ_sgl = single(cZ);
assert(isa(cZ_sgl.c,'single'));
assert(isa(cZ_sgl.G,'single'));

% test zonoBundle
Z1 = zonotope([1;0],[1 0; 0 1]);
Z2 = zonotope([0;0],[2 0; 0 2]);
zB = zonoBundle({Z1,Z2});
zB_sgl = single(zB);
assert(isa(zB_sgl.Z{1}.c,'single'));
assert(isa(zB_sgl.Z{2}.c,'single'));

% test emptySet
O = emptySet(2);
O_sgl = single(O);
assert(isa(O_sgl,'emptySet'));

% test fullspace
fs = fullspace(2);
fs_sgl = single(fs);
assert(isa(fs_sgl,'fullspace'));

% test that single of single keeps single
Z = zonotope(single([1;0]),single([1 0; 0 1]));
Z_sgl = single(Z);
assert(isa(Z_sgl.c,'single'));
assert(isa(Z_sgl.G,'single'));

% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
