function res = test_pointSegmentSignal_allTrue
% test_pointSegmentSignal_allTrue - unit test function of allTrue
%
% Syntax:
%    res = test_pointSegmentSignal_allTrue
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
sig = pointSegmentSignal([0,1,2,2.5,3,4], [tt,tt,ff,tt,tt,ff,ff,tt,ff,ff,ff,tt]);

% test case definition
test_cases = {
    % {int, expected}
    {stlInterval(0,1,true,false), true};
    {stlInterval(0,1,true,true), false};
    {stlInterval(0,2,true,false), false};
    {stlInterval(2,2.25,false,true), false};
    {stlInterval(0,inf), false};
    {stlInterval(4,inf,false,false), true};
    {stlInterval(100,inf), true};
};

% run tests
for i = 1:length(test_cases)
    int = test_cases{i}{1};
    expected = test_cases{i}{2};
    actual = sig.allTrue(int);
    assertLoop(actual == expected,i)
end

test_cases_neg = {
    % {int, expected}
    {stlInterval(3,4,true,true), true};
    {stlInterval(1,2.25,true,false), false};
};

% check negated signal
for i = 1:length(test_cases_neg)
    int = test_cases_neg{i}{1};
    expected = test_cases_neg{i}{2};
    actual = sig.allTrue(int,@(x) ~x);
    assertLoop(actual == expected,i)
end

res = true;

% ------------------------------ END OF CODE ------------------------------
