function res = test_kleeneSignal_or
% test_kleeneSignal_or - unit test function of or
%
% Syntax:
%    res = test_kleeneSignal_or
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
        switch kleeneVals(i) | kleeneVals(j)
            case kleene.True
                test_cases{end}.expected = {
                    interval;
                    [interval.toLeft(),interval.toRight()];
                    repmat(stlInterval(),0);
                };
            case kleene.Unknown
                test_cases{end}.expected = {
                    repmat(stlInterval(),0);
                    stlInterval(0,inf);
                    repmat(stlInterval(),0);
                };
            case kleene.False
                test_cases{end}.expected = {
                    repmat(stlInterval(),0);
                    [interval.toLeft(),interval.toRight()];
                    interval;
                };
        end
    end
end

test_cases{end+1}.lhs = kleeneSignal.indicator(stlInterval(2,4),tt,ff);
test_cases{end}.rhs = kleeneSignal.indicator(stlInterval(3,5,false,true),uu,ff);
trueInt = stlInterval(2,4);
unkInt = stlInterval(2,4).toRight() & stlInterval(3,5,false,true);
test_cases{end}.expected = {
    trueInt;
    unkInt;
    [trueInt.toLeft(),unkInt.toRight()];
};

for i = 1:length(test_cases)
    lhs = test_cases{i}.lhs;
    rhs = test_cases{i}.rhs;
    expected = test_cases{i}.expected;
    actual = lhs | rhs;
    actualInt = arrayfun(@(match) actual.findIntervals(match),kleeneVals,'UniformOutput',false);

    assertLoop(compareIntervals(actualInt{1},expected{1}),i)
    assertLoop(compareIntervals(actualInt{2},expected{2}),i)
    assertLoop(compareIntervals(actualInt{3},expected{3}),i)
end

res = true;

% ------------------------------ END OF CODE ------------------------------
