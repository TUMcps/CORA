function res = test_fullspace_radius
% test_fullspace_radius - unit test function of radius
%
% Syntax:
%    res = test_fullspace_radius
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
val = radius(fs);

% check result
res = val == Inf;

% ------------------------------ END OF CODE ------------------------------
