function res = test_spectraShadow_isemptyobject
% test_spectraShadow_isemptyobject - unit test function of isemptyobject
%
% Syntax:
%    res = test_spectraShadow_isemptyobject
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

% Authors:       Adrian Kulmburg
% Written:       14-August-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------


% empty spectrahedron
SpS_empty = spectraShadow.empty();
assert(isemptyobject(SpS_empty));


% 1D, bounded, non-degenerate
SpS = spectraShadow([1 0 1 0;0 1 0 -1]);
assert(~isemptyobject(SpS));

% 1D, empty
SpS = spectraShadow([-1 0]);
assert(~isemptyobject(SpS));

% 1D, unbounded
SpS = spectraShadow([1 0]);
assert(~isemptyobject(SpS));

% 1D, single point
SpS = spectraShadow([-1 0 1 0;0 1 0 -1]);
assert(~isemptyobject(SpS));


% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
