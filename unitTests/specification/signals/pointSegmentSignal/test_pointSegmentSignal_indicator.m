function res = test_pointSegmentSignal_indicator
% test_pointSegmentSignal_indicator - unit test function of indicator
%
% Syntax:
%    res = test_pointSegmentSignal_indicator
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

tt = true;
ff = false;

% test case definition
test_cases = {
    % {int, val, default, expected}
    {stlInterval(3,4,true,true),tt,ff,pointSegmentSignal([0,3,4],[ff,ff,tt,tt,tt,ff])};
    {stlInterval(3,4,false,true),tt,ff,pointSegmentSignal([0,3,4],[ff,ff,ff,tt,tt,ff])};
    {stlInterval(3,4,true,false),tt,ff,pointSegmentSignal([0,3,4],[ff,ff,tt,tt,ff,ff])};
    {stlInterval(3,4,false,false),tt,ff,pointSegmentSignal([0,3,4],[ff,ff,ff,tt,ff,ff])};
    {stlInterval(3,4,true,true),ff,tt,pointSegmentSignal([0,3,4],[tt,tt,ff,ff,ff,tt])};
    {stlInterval(3,4,true,true),tt,tt,pointSegmentSignal(0,[tt,tt])};
    {stlInterval(3,4,true,true),ff,ff,pointSegmentSignal(0,[ff,ff])};
    {stlInterval(),tt,tt,pointSegmentSignal(0,[tt,tt])};
    {stlInterval(),tt,ff,pointSegmentSignal(0,[ff,ff])};
    {stlInterval(0,inf),tt,tt,pointSegmentSignal(0,[tt,tt])};
    {stlInterval(0,inf),tt,ff,pointSegmentSignal(0,[tt,tt])};
};

% run tests
for i = 1:length(test_cases)
    int = test_cases{i}{1};
    val = test_cases{i}{2};
    default = test_cases{i}{3};
    expected = test_cases{i}{4};
    actual = pointSegmentSignal.indicator(int,val,default);
    assertLoop(actual == expected,i)
end

res = true;

% ------------------------------ END OF CODE ------------------------------
