function res = test_fourValuedSignal_combine
% test_fourValuedSignal_combine - unit test function of combine
%
% Syntax:
%    res = test_fourValuedSignal_combine
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
        switch (~fourValuedVals(i) | fourValuedVals(j))
            case fourValued.True
                % case 1
                test_cases{end}.expected = {
                    interval;
                    repmat(stlInterval(),0);
                    repmat(stlInterval(),0);
                    [interval.toLeft(),interval.toRight()];
                };
            case fourValued.Unknown
                % case 2
                test_cases{end}.expected = {
                    repmat(stlInterval(),0);
                    interval;
                    repmat(stlInterval(),0);
                    [interval.toLeft(),interval.toRight()];
                };
            case fourValued.False
                % case 3
                test_cases{end}.expected = {
                    repmat(stlInterval(),0);
                    repmat(stlInterval(),0);
                    interval;
                    [interval.toLeft(),interval.toRight()];
                };
            case fourValued.Inconclusive
                % case 4
                test_cases{end}.expected = {
                    repmat(stlInterval(),0);
                    repmat(stlInterval(),0);
                    repmat(stlInterval(),0);
                    stlInterval(0,inf);
                };
        end
    end
end

% case 5
test_cases{end+1}.lhs = fourValuedSignal.indicator(stlInterval(2,4),ff,tt);
test_cases{end}.rhs = fourValuedSignal.indicator(stlInterval(3,5,false,true),ii,uu);
trueInt = stlInterval(2,4);
incInt = stlInterval(2,4).toRight() & stlInterval(3,5,false,true);
test_cases{end}.expected = {
    trueInt;
    [trueInt.toLeft(),incInt.toRight()];
    repmat(stlInterval(),0);
    incInt;
};

% test cases
for i = 1:length(test_cases)
    lhs = test_cases{i}.lhs;
    rhs = test_cases{i}.rhs;
    expected = test_cases{i}.expected;
    actual = fourValuedSignal.combine(@aux_implies_and_neg_lhs,lhs,rhs);
    impliesInt = arrayfun(@(match) actual(1).findIntervals(match),fourValuedVals,'UniformOutput',false);
    negInt = arrayfun(@(match) actual(2).findIntervals(match),fourValuedVals,'UniformOutput',false);
    negIntExpected = arrayfun(@(match) findIntervals(~lhs,match),fourValuedVals,'UniformOutput',false);
    
    % check whether we got the expected values
    assertLoop(compareIntervals(impliesInt{1},expected{1}),i)
    assertLoop(compareIntervals(impliesInt{2},expected{2}),i)
    assertLoop(compareIntervals(impliesInt{3},expected{3}),i)
    assertLoop(compareIntervals(impliesInt{4},expected{4}),i)
    assertLoop(compareIntervals(negInt{1},negIntExpected{1}),i)
    assertLoop(compareIntervals(negInt{2},negIntExpected{2}),i)
    assertLoop(compareIntervals(negInt{3},negIntExpected{3}),i)
    assertLoop(compareIntervals(negInt{4},negIntExpected{4}),i)
end

res = true;
end


% Auxiliary functions -----------------------------------------------------

function out = aux_implies_and_neg_lhs(args)
    out(1) = ~args(1) | args(2);
    out(2) = ~args(1);
end

% ------------------------------ END OF CODE ------------------------------
