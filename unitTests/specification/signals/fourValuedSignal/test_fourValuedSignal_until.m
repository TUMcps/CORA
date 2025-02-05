function res = test_fourValuedSignal_until
% test_fourValuedSignal_until - unit test function of until
%
% Syntax:
%    res = test_fourValuedSignal_until
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
fourValuedVals = [tt,uu,ff,ii];
test_cases = {};

% case 1
test_cases{end+1}.lhs = fourValuedSignal.indicator(stlInterval(0,5,true,false),tt,ii);
test_cases{end}.interval = stlInterval(0,1);
test_cases{end}.rhs = fourValuedSignal.indicator(stlInterval(0,1,true,false),uu,ii) ...
    .set(stlInterval(1,1.5,true,false),tt) ...
    .set(stlInterval(1.5,2,true,false),ff) ...
    .set(stlInterval(2,3,true,false),tt) ...
    .set(stlInterval(3,5,true,false),ff);
test_cases{end}.expected = {
    [stlInterval(0,3,true,false)];
    repmat(stlInterval(),0);
    [stlInterval(3,4,true,false)];
    [stlInterval(4,inf)];
};

% case 2
test_cases{end+1}.lhs = fourValuedSignal.indicator(stlInterval(0,1,true,false),uu,ii) ...
    .set(stlInterval(1,5,true,false),tt);
test_cases{end}.interval = stlInterval(0,1);
test_cases{end}.rhs = fourValuedSignal.indicator(stlInterval(0,1,true,false),uu,ii) ...
    .set(stlInterval(1,1.5,true,false),tt) ...
    .set(stlInterval(1.5,2,true,false),ff) ...
    .set(stlInterval(2,3,true,false),tt) ...
    .set(stlInterval(3,5,true,false),ff);
test_cases{end}.expected = {
    [stlInterval(1,3,true,false)];
    [stlInterval(0,1,true,false)];
    [stlInterval(3,4,true,false)];
    [stlInterval(4,inf)];
};

% case 3
test_cases{end+1}.lhs = fourValuedSignal.indicator(stlInterval(0,4,true,false),tt,ii) ...
    .set(stlInterval(4,5,true,false),uu);
test_cases{end}.interval = stlInterval(0,1);
test_cases{end}.rhs = fourValuedSignal.indicator(stlInterval(0,1,true,false),uu,ii) ...
    .set(stlInterval(1,1.5,true,false),tt) ...
    .set(stlInterval(1.5,2,true,false),ff) ...
    .set(stlInterval(2,3,true,false),tt) ...
    .set(stlInterval(3,5,true,false),ff);
test_cases{end}.expected = {
    [stlInterval(0,3,true,false)];
    repmat(stlInterval(),0);
    [stlInterval(3,4,true,false)];
    [stlInterval(4,inf)];
};

% case 4
test_cases{end+1}.lhs = fourValuedSignal.indicator(stlInterval(0,4,true,false),tt,ii) ...
    .set(stlInterval(4,5,true,false),ff);
test_cases{end}.interval = stlInterval(0,1);
test_cases{end}.rhs = fourValuedSignal.indicator(stlInterval(0,1,true,false),uu,ii) ...
    .set(stlInterval(1,1.5,true,false),tt) ...
    .set(stlInterval(1.5,2,true,false),ff) ...
    .set(stlInterval(2,3,true,false),tt) ...
    .set(stlInterval(3,5,true,false),ff);
test_cases{end}.expected = {
    [stlInterval(0,3,true,false)];
    repmat(stlInterval(),0);
    [stlInterval(3,5,true,false)];
    [stlInterval(5,inf)];
};

% case 5
test_cases{end+1}.lhs = fourValuedSignal.indicator(stlInterval(0,3,false),tt,ii) ...
    .set(stlInterval(3,5,true,false),uu) ...
    .set(stlInterval(5,6,true,false),ff);
test_cases{end}.interval = stlInterval(0,1,false);
test_cases{end}.rhs = fourValuedSignal.indicator(stlInterval(0,2,true,false),uu,ii) ...
    .set(stlInterval(2,4,true,false),tt);
test_cases{end}.expected = {
    [stlInterval(1,3,false)];
    [stlInterval(0,1,true),stlInterval(3,4,true,false)];
    [stlInterval(5,6,true,false)];
    [stlInterval(4,5,true,false),stlInterval(6,inf)];
    % inconclusive in [4, 5), since rhs is inconclusive in (4, 6)
    % thus it could still become false in this interval,
    % which would make until false in [4, 5)
};

% case 6
test_cases{end+1}.lhs = fourValuedSignal.indicator(stlInterval(0,3,false),tt,ii) ...
    .set(stlInterval(3,5,true,false),uu) ...
    .set(stlInterval(5,6,true,false),ff);
test_cases{end}.interval = stlInterval(1,1.5);
test_cases{end}.rhs = fourValuedSignal.indicator(stlInterval(0,2,true,false),uu,ii) ...
    .set(stlInterval(2,4,true,false),tt);
test_cases{end}.expected = {
    [stlInterval(0.5,2,true)];
    [stlInterval(0,0.5,true,false),stlInterval(2,3,false)];
    [stlInterval(4,6,false)];
    [stlInterval(3,4,true),stlInterval(6,inf)];
    % inconclusive in [3, 4], since rhs is inconclusive in [4, 5.5]
    % thus it could still become false in this interval,
    % which would make until false in [3, 4]
};

% case 7
test_cases{end+1}.lhs = fourValuedSignal.indicator(stlInterval(0,1,true,false),uu,ii) ...
    .set(stlInterval(1),ff);
test_cases{end}.interval = stlInterval(0,2,false);
test_cases{end}.rhs = fourValuedSignal.indicator(stlInterval(0,1),tt,ii);
test_cases{end}.expected = {
    repmat(stlInterval(),0);
    [stlInterval(0,1,true,false)];
    repmat(stlInterval(),0);
    [stlInterval(1,inf)];
    % unknown in [0, 1), since rhs is unknown in [0, 1) and false at 1
    % thus until becomes false once we cross time 1
    % however, since rhs is true in [0, 1], until is unknown if we do not cross 1
    % thus, we can say with certainty that until is unknown,
    % even though some parts of its future reach are still inconclusive
};

% test cases
for i = 1:length(test_cases)
    lhs = test_cases{i}.lhs;
    interval = test_cases{i}.interval;
    rhs = test_cases{i}.rhs;
    expected = test_cases{i}.expected;
    actual = until(lhs,interval,rhs);
    actualInt = arrayfun(@(match) actual.findIntervals(match),fourValuedVals,'UniformOutput',false);
    
    assertLoop(compareIntervals(actualInt{1},expected{1}),i)
    assertLoop(compareIntervals(actualInt{2},expected{2}),i)
    assertLoop(compareIntervals(actualInt{3},expected{3}),i)
    assertLoop(compareIntervals(actualInt{4},expected{4}),i)
end

res = true;

% ------------------------------ END OF CODE ------------------------------
