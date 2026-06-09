function res = test_contSet_sparse
% test_contSet_sparse - unit test function of sparse
%
% Syntax:
%    res = test_contSet_sparse
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
I_sp = sparse(I);
assert(issparse(I_sp.inf));
assert(issparse(I_sp.sup));

% test zonotope
Z = zonotope([1;0],[1 0; 0 1]);
Z_sp = sparse(Z);
assert(issparse(Z_sp.c));
assert(issparse(Z_sp.G));

% test capsule
C = capsule([1;0],[1;1],0.5);
C_sp = sparse(C);
assert(issparse(C_sp.c));
assert(issparse(C_sp.g));
assert(issparse(C_sp.r));

% test ellipsoid
E = ellipsoid([1 0; 0 2],[1;0]);
E_sp = sparse(E);
assert(issparse(E_sp.Q));
assert(issparse(E_sp.q));

% test polytope (H-rep)
P = polytope([1 0; -1 0; 0 1; 0 -1],[1;1;1;1]);
P_sp = sparse(P);
assert(issparse(P_sp.A));
assert(issparse(P_sp.b));

% test polytope (V-rep)
P = polytope([1 0; 0 1; -1 0]);
P_sp = sparse(P);
assert(issparse(P_sp.V));

% test polyZonotope
pZ = polyZonotope([1;0],[1 0; 0 1],[0.5;0],[1 0; 0 1]);
pZ_sp = sparse(pZ);
assert(issparse(pZ_sp.c));
assert(issparse(pZ_sp.G));
assert(issparse(pZ_sp.GI));

% test conZonotope
cZ = conZonotope([0;0],[1 0; 0 1],[1 -1],1);
cZ_sp = sparse(cZ);
assert(issparse(cZ_sp.c));
assert(issparse(cZ_sp.G));

% test zonoBundle
Z1 = zonotope([1;0],[1 0; 0 1]);
Z2 = zonotope([0;0],[2 0; 0 2]);
zB = zonoBundle({Z1,Z2});
zB_sp = sparse(zB);
assert(issparse(zB_sp.Z{1}.c));
assert(issparse(zB_sp.Z{2}.c));

% test emptySet
O = emptySet(2);
O_sp = sparse(O);
assert(isa(O_sp,'emptySet'));

% test fullspace
fs = fullspace(2);
fs_sp = sparse(fs);
assert(isa(fs_sp,'fullspace'));

% test that sparse of sparse keeps sparse
Z = zonotope(sparse([1;0]),sparse([1 0; 0 1]));
Z_sp = sparse(Z);
assert(issparse(Z_sp.c));
assert(issparse(Z_sp.G));

% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
