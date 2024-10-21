function res = test_fourValuedSignal_and
% test_fourValuedSignal_and - unit test function of and
%
% Syntax:
%    res = test_fourValuedSignal_and
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
% Written:       20-February-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% signals
tt = fourValued.True;
uu = fourValued.Unknown;
ff = fourValued.False;
ii = fourValued.Inconclusive;

interval = stlInterval(1,2);
fourValuedVals = [tt,uu,ff,ii];
test_cases = {};
for i = 1:length(fourValuedVals)
    for j = 1:length(fourValuedVals)
        test_cases{end+1}.lhs = fourValuedSignal.indicator(interval,fourValuedVals(i),ii);
        test_cases{end}.rhs = fourValuedSignal.indicator(interval,fourValuedVals(j),ii);
        switch fourValuedVals(i) & fourValuedVals(j)
            case fourValued.True
                test_cases{end}.expected = {
                    interval;
                    repmat(stlInterval(),0);
                    repmat(stlInterval(),0);
                    [interval.toLeft(),interval.toRight()];
                };
            case fourValued.Unknown
                test_cases{end}.expected = {
                    repmat(stlInterval(),0);
                    interval;
                    repmat(stlInterval(),0);
                    [interval.toLeft(),interval.toRight()];
                };
            case fourValued.False
                test_cases{end}.expected = {
                    repmat(stlInterval(),0);
                    repmat(stlInterval(),0);
                    interval;
                    [interval.toLeft(),interval.toRight()];
                };
            case fourValued.Inconclusive
                test_cases{end}.expected = {
                    repmat(stlInterval(),0);
                    repmat(stlInterval(),0);
                    repmat(stlInterval(),0);
                    stlInterval(0,inf);
                };
        end
    end
end

test_cases{end+1}.lhs = fourValuedSignal.indicator(stlInterval(2,4),ff,tt);
test_cases{end}.rhs = fourValuedSignal.indicator(stlInterval(3,5,false,true),ii,uu);
falseInt = stlInterval(2,4);
incInt = stlInterval(2,4).toRight() & stlInterval(3,5,false,true);
test_cases{end}.expected = {
    repmat(stlInterval(),0);
    [falseInt.toLeft(),incInt.toRight()];
    falseInt;
    incInt;
};

for i = 1:length(test_cases)
    lhs = test_cases{i}.lhs;
    rhs = test_cases{i}.rhs;
    expected = test_cases{i}.expected;
    actual = lhs & rhs;
    actualInt = arrayfun(@(match) actual.findIntervals(match),fourValuedVals,'UniformOutput',false);
    
    assertLoop(compareIntervals(actualInt{1},expected{1}),i)
    assertLoop(compareIntervals(actualInt{2},expected{2}),i)
    assertLoop(compareIntervals(actualInt{3},expected{3}),i)
    assertLoop(compareIntervals(actualInt{4},expected{4}),i)
end

res = true;

% ------------------------------ END OF CODE ------------------------------
