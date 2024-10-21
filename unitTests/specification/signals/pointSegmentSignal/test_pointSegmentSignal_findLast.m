function res = test_pointSegmentSignal_findLast
% test_pointSegmentSignal_findLast - unit test function of findLast
%
% Syntax:
%    res = test_pointSegmentSignal_findLast
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

% Authors:       Florian Lercher
% Written:       16-February-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% signals
tt = true;
ff = false;

% test case definition
test_cases = {
    % {sig, expTimeT, expInclT, expTimeF, expInclF}
    {pointSegmentSignal(0,[tt,ff]), 0, true, inf, false};
    {pointSegmentSignal(0,[tt,tt]), inf, false, 0, false};
    {pointSegmentSignal([0,1],[tt,ff,tt,ff]), 1, true, inf, false};
    {pointSegmentSignal([0,1],[tt,tt,ff,ff]), 1, false, inf, false};
};

% run tests
for i = 1:length(test_cases)
    sig = test_cases{i}{1};
    expTimeT = test_cases{i}{2};
    expInclT = test_cases{i}{3};
    expTimeF = test_cases{i}{4};
    expInclF = test_cases{i}{5};
    [actTimeT,actInclT] = sig.findLast(true);
    [actTimeF,actInclF] = sig.findLast(false);
    
    assertLoop(actTimeT == expTimeT,i)
    assertLoop(actInclT == expInclT,i)
    assertLoop(actTimeF == expTimeF,i)
    assertLoop(actInclF == expInclF,i)
end

res = true;

% ------------------------------ END OF CODE ------------------------------
