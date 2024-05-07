function res = test_polytope_mpt
% test_polytope_mpt - unit test function for various polytopes functions,
%    where the evaluation using the mpt implementation returns wrong
%    results
%
% Syntax:
%    res = test_polytope_mpt
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

% Authors:       Mark Wetzlinger
% Written:       02-May-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true(0);

% vertices
A = [1 0; -1 1; -1 -1];
b = [1;1;1];
P = polytope(A,b);
% P_mpt = Polyhedron(A,b);
% P_mpt.V % wrong when called without cddmex on the path

V = vertices(P);
V_true = [1 2; -1 0; 1 -2]';
res(end+1,1) = compareMatrices(V,V_true);

% center
% Any unbounded polytope returns a  large value as center, instead of NAN
P_unb = polytope([1 -0.1; 0.1 -1; -0.1 -1; -1 -0.1],ones(4,1));
c = center(P);
% Ph_unb = Polyhedron([1 -0.1; 0.1 -1; -0.1 -1; -1 -0.1],ones(4,1));
% ch= chebyCenter(Ph_unb).x

% Eq/vertices
% If you compute the vertices of an unbounded polytope, it starts returning false results for equality
A = [0 1 0; 0 0 1; 0 -1 0; 0 0 -1]; b = ones(4,1); Ae = [1 0 0]; be = 2;
P = polytope(A,b,Ae,be);
%
% Ph_d_unb = Polyhedron('A', A,'b', b, 'Ae', Ae, 'be', be);
% Ph_d_unb == Ph_d_unb % will be correct
% Ph_d_unb.V
% Ph_d_unb == Ph_d_unb % now result will be wrong


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
