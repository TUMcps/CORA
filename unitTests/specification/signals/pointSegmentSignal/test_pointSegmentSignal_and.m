function res = test_pointSegmentSignal_and
% test_pointSegmentSignal_and - unit test function of and
%
% Syntax:
%    res = test_pointSegmentSignal_and
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

test_cases{1}.lhs = pointSegmentSignal([0,1],[tt,ff,tt,ff]);
test_cases{1}.rhs = pointSegmentSignal([0,1],[tt,tt,ff,ff]);
test_cases{1}.exp = pointSegmentSignal(0,[tt,ff]);

test_cases{2}.lhs = pointSegmentSignal([0,5],[tt,tt,ff,tt]);
test_cases{2}.rhs = pointSegmentSignal([0,2],[tt,tt,ff,tt]);
test_cases{2}.exp = pointSegmentSignal([0,2,5],[tt,tt,ff,tt,ff,tt]);

test_cases{3}.lhs = pointSegmentSignal([0,5],[ff,ff,tt,ff]);
test_cases{3}.rhs = pointSegmentSignal([0,2],[tt,tt,ff,tt]);
test_cases{3}.exp = pointSegmentSignal([0,5],[ff,ff,tt,ff]);

test_cases{4}.lhs = pointSegmentSignal([0,5],[ff,ff,tt,ff]);
test_cases{4}.rhs = pointSegmentSignal([0,2],[ff,ff,tt,ff]);
test_cases{4}.exp = pointSegmentSignal(0,[ff,ff]);

% run tests
for i = 1:length(test_cases)
    c = test_cases{i};
    expected = c.exp;
    actual = c.lhs & c.rhs;
    assertLoop(actual == expected,i)
end

% check multiple signals
for i = 1:length(test_cases)
    c = test_cases{i};
    expected = c.exp;
    actual = pointSegmentSignal.and_(c.lhs,c.rhs,c.exp);
    assertLoop(actual == expected,i)
end

res = true;

% ------------------------------ END OF CODE ------------------------------
