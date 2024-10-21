function res = test_kleeneSignal_until
% test_kleeneSignal_until - unit test function of until
%
% Syntax:
%    res = test_kleeneSignal_until
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
tt = kleene.True;
uu = kleene.Unknown;
ff = kleene.False;
kleeneVals = [tt,uu,ff];
test_cases = {};

test_cases{end+1}.lhs = kleeneSignal.indicator(stlInterval(0,5,true,false),tt,ff);
test_cases{end}.interval = stlInterval(0,1);
test_cases{end}.rhs = kleeneSignal.indicator(stlInterval(0,1,true,false),uu,ff) ...
    .set(stlInterval(1,1.5,true,false),tt) ...
    .set(stlInterval(1.5,2,true,false),ff) ...
    .set(stlInterval(2,3,true,false),tt) ...
    .set(stlInterval(3,5,true,false),ff);
test_cases{end}.expected = {
    [stlInterval(0,3,true,false)];
    repmat(stlInterval(),0);
    [stlInterval(3,inf)];
};

test_cases{end+1}.lhs = kleeneSignal.indicator(stlInterval(0,1,true,false),uu,ff) ...
    .set(stlInterval(1,5,true,false),tt);
test_cases{end}.interval = stlInterval(0,1);
test_cases{end}.rhs = kleeneSignal.indicator(stlInterval(0,1,true,false),uu,ff) ...
    .set(stlInterval(1,1.5,true,false),tt) ...
    .set(stlInterval(1.5,2,true,false),ff) ...
    .set(stlInterval(2,3,true,false),tt) ...
    .set(stlInterval(3,5,true,false),ff);
test_cases{end}.expected = {
    [stlInterval(1,3,true,false)];
    [stlInterval(0,1,true,false)];
    [stlInterval(3,inf)];
};

test_cases{end+1}.lhs = kleeneSignal.indicator(stlInterval(0,4,true,false),tt,ff) ...
    .set(stlInterval(4,5,true,false),uu);
test_cases{end}.interval = stlInterval(0,1);
test_cases{end}.rhs = kleeneSignal.indicator(stlInterval(0,1,true,false),uu,ff) ...
    .set(stlInterval(1,1.5,true,false),tt) ...
    .set(stlInterval(1.5,2,true,false),ff) ...
    .set(stlInterval(2,3,true,false),tt) ...
    .set(stlInterval(3,5,true,false),ff);
test_cases{end}.expected = {
    [stlInterval(0,3,true,false)];
    repmat(stlInterval(),0);
    [stlInterval(3,inf)];
};

test_cases{end+1}.lhs = kleeneSignal.indicator(stlInterval(0,4,true,false),tt,ff);
test_cases{end}.interval = stlInterval(0,1);
test_cases{end}.rhs = kleeneSignal.indicator(stlInterval(0,1,true,false),uu,ff) ...
    .set(stlInterval(1,1.5,true,false),tt) ...
    .set(stlInterval(1.5,2,true,false),ff) ...
    .set(stlInterval(2,3,true,false),tt) ...
    .set(stlInterval(3,5,true,false),ff);
test_cases{end}.expected = {
    [stlInterval(0,3,true,false)];
    repmat(stlInterval(),0);
    [stlInterval(3,inf)];
};

test_cases{end+1}.lhs = kleeneSignal.indicator(stlInterval(0,3,false),tt,ff) ...
    .set(stlInterval(3,5,true,false),uu);
test_cases{end}.interval = stlInterval(0,1,false);
test_cases{end}.rhs = kleeneSignal.indicator(stlInterval(0,2,true,false),uu,ff) ...
    .set(stlInterval(2,4,true,false),tt);
test_cases{end}.expected = {
    [stlInterval(1,3,false)];
    [stlInterval(0,1,true),stlInterval(3,4,true,false)];
    [stlInterval(4,inf)];
};

test_cases{end+1}.lhs = kleeneSignal.indicator(stlInterval(0,3,false),tt,ff) ...
    .set(stlInterval(3,5,true,false),uu);
test_cases{end}.interval = stlInterval(1,1.5);
test_cases{end}.rhs = kleeneSignal.indicator(stlInterval(0,2,true,false),uu,ff) ...
    .set(stlInterval(2,4,true,false),tt);
test_cases{end}.expected = {
    [stlInterval(0.5,2,true)];
    [stlInterval(0,0.5,true,false),stlInterval(2,3,false)];
    [stlInterval(3,inf)];
};

for i = 1:length(test_cases)
    lhs = test_cases{i}.lhs;
    interval = test_cases{i}.interval;
    rhs = test_cases{i}.rhs;
    expected = test_cases{i}.expected;
    actual = until(lhs,interval,rhs);
    actualInt = arrayfun(@(match) actual.findIntervals(match),kleeneVals,'UniformOutput',false);
    
    assertLoop(compareIntervals(actualInt{1},expected{1}),i)
    assertLoop(compareIntervals(actualInt{2},expected{2}),i)
    assertLoop(compareIntervals(actualInt{3},expected{3}),i)
end

res = true;

% ------------------------------ END OF CODE ------------------------------
