function res = test_pointSegmentSignal_at
% test_pointSegmentSignal_at - unit test function of at
%
% Syntax:
%    res = test_pointSegmentSignal_at
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
sig = pointSegmentSignal([0,1,2,2.5], [tt,tt,ff,tt,tt,ff,ff,tt]);

% test case definition
test_cases = {
    % {time, expected}
    {0, tt};
    {0.5, tt};
    {1, ff};
    {1.5, tt};
    {2, tt};
    {2.25, ff};
    {2.5, ff};
    {3, tt};
    {inf, tt};
};

% run tests
for i = 1:length(test_cases)
    time = test_cases{i}{1};
    expected = test_cases{i}{2};
    actual = sig.at(time);
    assertLoop(actual == expected,i)
end

res = true;

% ------------------------------ END OF CODE ------------------------------
