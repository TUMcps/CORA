function res = test_emptySet_volume
% test_emptySet_volume - unit test function of volume
%
% Syntax:
%    res = test_emptySet_volume
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

% compute volume
val = volume(O);
res = val == 0;

% ------------------------------ END OF CODE ------------------------------
