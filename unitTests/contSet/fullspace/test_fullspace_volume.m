function res = test_fullspace_volume
% test_fullspace_volume - unit test function of volume
%
% Syntax:
%    res = test_fullspace_volume
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

% init fullspace
n = 2;
fs = fullspace(n);

% compute volume
val = volume(fs);

% check result
res = val == Inf;

% ------------------------------ END OF CODE ------------------------------
