function res = test_fullspace_isFullDim
% test_fullspace_isFullDim - unit test function of isFullDim
%
% Syntax:
%    res = test_fullspace_isFullDim
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

% check property
res = isFullDim(fs);

% ------------------------------ END OF CODE ------------------------------
