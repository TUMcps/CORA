function res = test_emptySet_center
% test_emptySet_center - unit test function of center
%
% Syntax:
%    res = test_emptySet_center
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
% See also: -

% Authors:       Mark Wetzlinger
% Written:       05-April-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% init empty set
n = 2;
O = emptySet(n);

% center
c = center(O);
assert(all(size(c) == [n,0]));

% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
