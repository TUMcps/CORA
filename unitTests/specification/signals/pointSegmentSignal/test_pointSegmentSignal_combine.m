function res = test_pointSegmentSignal_combine
% test_pointSegmentSignal_combine - unit test function of combine
%
% Syntax:
%    res = test_pointSegmentSignal_combine
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

test_cases{1}.lhs = pointSegmentSignal([0,1],[ff,tt,ff,tt]);
test_cases{1}.rhs = pointSegmentSignal([0,1],[tt,tt,ff,ff]);
test_cases{1}.exp = pointSegmentSignal([0,1],[tt,tt,tt,ff]);

test_cases{2}.lhs = pointSegmentSignal([0,5],[ff,ff,tt,ff]);
test_cases{2}.rhs = pointSegmentSignal([0,2],[tt,tt,ff,tt]);
test_cases{2}.exp = pointSegmentSignal(0,[tt,tt]);

test_cases{3}.lhs = pointSegmentSignal([0,5],[tt,tt,ff,tt]);
test_cases{3}.rhs = pointSegmentSignal([0,2],[tt,tt,ff,tt]);
test_cases{3}.exp = pointSegmentSignal([0,2],[tt,tt,ff,tt]);

test_cases{4}.lhs = pointSegmentSignal([0,5],[tt,tt,ff,tt]);
test_cases{4}.rhs = pointSegmentSignal([0,2],[ff,ff,tt,ff]);
test_cases{4}.exp = pointSegmentSignal([0,2,5],[ff,ff,tt,ff,tt,ff]);

% run tests
for i = 1:length(test_cases)
    c = test_cases{i};
    expected = c.exp;
    actual = pointSegmentSignal.combine(@aux_implies_and_neg_lhs,c.lhs,c.rhs);
    assertLoop(actual(1) == expected,i)
    assertLoop(actual(2) == ~c.lhs,i)
end

res = true;
end


% Auxiliary functions -----------------------------------------------------

function out = aux_implies_and_neg_lhs(args)
    out(1) = ~args(1) | args(2);
    out(2) = ~args(1);
end

% ------------------------------ END OF CODE ------------------------------
