function res = test_pointSegmentSignal_cutAtFirstFallingEdge
% test_pointSegmentSignal_cutAtFirstFallingEdge - unit test function of cutAtFirstFallingEdge
%
% Syntax:
%    res = test_pointSegmentSignal_cutAtFirstFallingEdge
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
    % {sig, expected}
    {pointSegmentSignal(0,[tt,ff]), pointSegmentSignal(0,[tt,ff])};
    {pointSegmentSignal(0,[ff,tt]), pointSegmentSignal(0,[ff,tt])}; % no falling edge here
    {pointSegmentSignal([0,1,2],[ff,tt,tt,ff,tt,tt]), pointSegmentSignal([0,1],[ff,tt,tt,ff])};
    {pointSegmentSignal([0,1,2],[ff,tt,ff,tt,ff,tt]), pointSegmentSignal([0,1],[ff,tt,ff,ff])};
    {pointSegmentSignal([0,2,3],[tt,tt,ff,ff,ff,tt]), pointSegmentSignal([0,2],[tt,tt,ff,ff])}
};

% run tests
for i = 1:length(test_cases)
    sig = test_cases{i}{1};
    expected = test_cases{i}{2};
    actual = sig.cutAtFirstFallingEdge();
    assertLoop(actual == expected,i)
end

res = true;

% ------------------------------ END OF CODE ------------------------------
