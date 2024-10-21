function res = test_spectraShadow_Inf
% test_spectraShadow_Inf - unit test function of R^n instantiation
%
% Syntax:
%    res = test_spectraShadow_Inf
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

% Authors:       Adrian Kulmburg
% Written:       13-June-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% 1D
n = 1;
SpS = spectraShadow.Inf(n);
assert(dim(SpS) == 1);

% 5D
n = 5;
SpS = spectraShadow.Inf(n);
assert(dim(SpS) == 5);

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
