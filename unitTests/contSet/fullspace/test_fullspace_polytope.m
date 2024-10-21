function res = test_fullspace_polytope
% test_fullspace_polytope - unit test function of polytope conversion
%
% Syntax:
%    res = test_fullspace_polytope
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
% Written:       15-December-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% init fullspace
fs = fullspace(2);
P = polytope(fs);
assert(~isBounded(P));
assert(supportFunc(P,[1;0],'upper') == Inf);
assert(supportFunc(P,[-1;0],'upper') == Inf);
assert(supportFunc(P,[0;1],'upper') == Inf);
assert(supportFunc(P,[0;-1],'upper') == Inf);


% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
