function res = test_spectraShadow_empty
% test_spectraShadow_empty - unit test function of empty instantiation
%
% Syntax:
%    res = test_spectraShadow_empty
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
SpS = spectraShadow.empty(n);
assert(representsa(SpS,'emptySet') && dim(SpS) == 1);

% 5D
n = 5;
SpS = spectraShadow.empty(n);
assert(representsa(SpS,'emptySet') && dim(SpS) == 5);

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
