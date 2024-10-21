function res = test_polytope_levelSet
% test_polytope_levelSet - unit test function of level set conversion
%
% Syntax:
%    res = test_polytope_levelSet
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
% Written:       15-December-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

syms x1 x2 x3;

% 2D, bounded
A = [1 0; -1 1; -1 -1]; b = [1;1;1];
P = polytope(A,b);
ls = levelSet(P);
eq1 = x1 - 1;
eq2 = -x1 + x2 - 1;
eq3 = -x1 - x2 - 1;
ls_true = levelSet([eq1;eq2;eq3],[x1;x2],{"<=","<=","<="});
assert(length(ls.eq) == 3 && length(ls.vars) == 2 && length(ls.compOp) == 3);
assert(isequal(ls,ls_true));

% 3D, unbounded
A = [1 1 -1]; b = -1;
P = polytope(A,b);
ls = levelSet(P);
eq1 = x1 + x2 - x3 + 1;
ls_true = levelSet(eq1,[x1;x2;x3],{"<="});
assert(length(ls.eq) == 1 && length(ls.vars) == 3 && length(ls.compOp) == 1);
assert(isequal(ls,ls_true));

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
