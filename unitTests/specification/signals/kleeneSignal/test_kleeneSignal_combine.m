function res = test_kleeneSignal_combine
% test_kleeneSignal_combine - unit test function of combine
%
% Syntax:
%    res = test_kleeneSignal_combine
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
% Written:       19-February-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% signals
tt = kleene.True;
uu = kleene.Unknown;
ff = kleene.False;

interval = stlInterval(1,2);
kleeneVals = [tt,uu,ff];
test_cases = {};
for i = 1:length(kleeneVals)
    for j = 1:length(kleeneVals)
        test_cases{end+1}.lhs = kleeneSignal.indicator(interval,kleeneVals(i),uu);
        test_cases{end}.rhs = kleeneSignal.indicator(interval,kleeneVals(j),uu);
        switch (~kleeneVals(i) | kleeneVals(j))
            case kleene.True
                % case 1
                test_cases{end}.expected = {
                    interval;
                    [interval.toLeft(),interval.toRight()];
                    repmat(stlInterval(),0);
                };
            case kleene.Unknown
                % case 2
                test_cases{end}.expected = {
                    repmat(stlInterval(),0);
                    stlInterval(0,inf);
                    repmat(stlInterval(),0);
                };
            case kleene.False
                % case 3
                test_cases{end}.expected = {
                    repmat(stlInterval(),0);
                    [interval.toLeft(),interval.toRight()];
                    interval;
                };
        end
    end
end

% case 4
test_cases{end+1}.lhs = kleeneSignal.indicator(stlInterval(2,4),ff,tt);
test_cases{end}.rhs = kleeneSignal.indicator(stlInterval(3,5,false,true),uu,ff);
trueInt = stlInterval(2,4);
unkInt = stlInterval(2,4).toRight() & stlInterval(3,5,false,true);
test_cases{end}.expected = {
    trueInt;
    unkInt;
    [trueInt.toLeft(),unkInt.toRight()];
};

% test cases
for i = 1:length(test_cases)
    lhs = test_cases{i}.lhs;
    rhs = test_cases{i}.rhs;
    expected = test_cases{i}.expected;
    actual = kleeneSignal.combine(@aux_implies_and_neg_lhs,lhs,rhs);
    impliesInt = arrayfun(@(match) actual(1).findIntervals(match),kleeneVals,'UniformOutput',false);
    negInt = arrayfun(@(match) actual(2).findIntervals(match),kleeneVals,'UniformOutput',false);
    negIntExpected = arrayfun(@(match) findIntervals(~lhs,match),kleeneVals,'UniformOutput',false);
    
    % check whether we got the expected values
    assertLoop(compareIntervals(impliesInt{1},expected{1}),i)
    assertLoop(compareIntervals(impliesInt{2},expected{2}),i)
    assertLoop(compareIntervals(impliesInt{3},expected{3}),i)
    assertLoop(compareIntervals(negInt{1},negIntExpected{1}),i)
    assertLoop(compareIntervals(negInt{2},negIntExpected{2}),i)
    assertLoop(compareIntervals(negInt{3},negIntExpected{3}),i)
end

res = true;
end


% Auxiliary functions -----------------------------------------------------

function out = aux_implies_and_neg_lhs(args)
    out(1) = ~args(1) | args(2);
    out(2) = ~args(1);
end

% ------------------------------ END OF CODE ------------------------------
