function res = testLong_spectraShadow_contains
% testLong_spectraShadow_contains - unit test function of contains
%
% Syntax:
%    res = testLong_spectraShadow_contains
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

runs = 5;

% Simple test: Spectrahedral shadow should contain its center-vector
for i=1:runs
    SpS = spectraShadow.generateRandom();
    c = center(SpS);
    assert(contains(SpS,c,'exact',1e-4))
end

res = true;

% The test of whether randPoint-points are contained in S will be done in
% testLong_spectraShadow_randPoint.m, so no need to do it twice

% ------------------------------ END OF CODE ------------------------------
