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

% setup
x = stl('x', 3);

% cases
p1 = finally(globally(x(1) < 5 & x(2) > 7, interval(1,2)), interval(2,3));
ap1 = atomicProposition(halfspace([1 0 0], 5));
ap2 = atomicProposition(halfspace([0 -1 0], -7));
r1 = finally(globally(stl('x1 < 5',ap1) & stl('x2 > 7',ap2), interval(1,2)), interval(2,3));

[e1,m1] = combineAtomicPropositions(p1);

res = isequal(r1, e1);
res(end+1,1) = isequal(ap1, m1("x1 < 5"));
res(end+1,1) = isequal(ap2, m1("x2 > 7"));

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
