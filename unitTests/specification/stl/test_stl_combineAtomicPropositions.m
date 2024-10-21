function res = test_stl_combineAtomicPropositions
% test_stl_combineAtomicPropositions - unit test of stl
%                                      combineAtomicPropositions method
%
% Syntax:
%    res = test_stl_combineAtomicPropositions
%
% Inputs:
%    -
%
% Outputs:
%    res - boolean
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Benedikt Seidl
% Written:       12-May-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;

% setup
x = stl('x', 3);

% cases
p1 = finally(globally(x(1) < 5 & x(2) > 7, interval(1,2)), interval(2,3));
ap1 = atomicProposition(polytope([1 0 0], 5));
ap2 = atomicProposition(polytope([0 -1 0], -7));
r1 = finally(globally(stl('x1 < 5',ap1) & stl('x2 > 7',ap2), interval(1,2)), interval(2,3));

[e1,m1] = combineAtomicPropositions(p1);

assert(isequal(r1, e1));
assert(isequal(ap1, m1("x1 < 5")));
assert(isequal(ap2, m1("x2 > 7")));

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
