function res = test_fullspace_copy
% test_fullspace_copy - unit test function of copy
%
% Syntax:
%    res = test_fullspace_copy
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

% 2D fullspace
fs = fullspace(2);
fs_copy = copy(fs);
assert(isequal(fs,fs_copy));


% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
