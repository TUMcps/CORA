function res = test_emptySet_copy
% test_emptySet_copy - unit test function of copy
%
% Syntax:
%    res = test_emptySet_copy
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

% 2D emptySet
O = emptySet(2);
O_copy = copy(O);
assert(isequal(O,O_copy));


% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
