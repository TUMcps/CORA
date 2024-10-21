function res = test_polytope_copy
% test_polytope_copy - unit test function of copy
%
% Syntax:
%    res = test_polytope_copy
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
% Written:       02-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% 2D polytope
P = polytope([1 0; -1 1; -1 -1],[1;1;1]);
P_copy = copy(P);

% ensure that set properties are independent from another
isFullDim(P);
assert(isempty(P_copy.fullDim.val));

% ensure that set properties are copied
P_copy = copy(P);
assert(P_copy.fullDim.val);

% ensure that sets are equal
assert(isequal(P,P_copy));

% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
