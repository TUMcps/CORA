function res = test_interval_copy
% test_interval_copy - unit test function of copy
%
% Syntax:
%    res = test_interval_copy
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

% 2D interval
I = interval(2,4);
I_copy = copy(I);
assert(isequal(I,I_copy));


% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
