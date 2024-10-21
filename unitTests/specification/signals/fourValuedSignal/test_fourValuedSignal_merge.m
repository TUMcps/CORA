function res = test_fourValuedSignal_merge
% test_fourValuedSignal_merge - unit test function of merge
%
% Syntax:
%    res = test_fourValuedSignal_merge
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

test_cases = {};

test_cases{end+1}.lhs = fourValuedSignal.indicator(stlInterval(0,2),tt,ii);
test_cases{end}.rhs = fourValuedSignal.indicator(stlInterval(3,5,false,true),uu,ii);
test_cases{end}.expected = test_cases{end}.lhs.set(stlInterval(3,5,false,true),uu);

test_cases{end+1}.lhs = fourValuedSignal.indicator(stlInterval(0,2),tt,ii);
test_cases{end}.rhs = fourValuedSignal.indicator(stlInterval(3,5,false,true),ii,uu);
test_cases{end}.expected = test_cases{end}.lhs.set(stlInterval(0,3),uu).set(stlInterval(5,inf,false),uu);

test_cases{end+1}.lhs = fourValuedSignal.indicator(stlInterval(0,2),tt,ii);
test_cases{end}.rhs = fourValuedSignal.indicator(stlInterval(3,5,false,true),tt,ff);
test_cases{end}.expected = test_cases{end}.rhs;

for i = 1:length(test_cases)
    lhs = test_cases{i}.lhs;
    rhs = test_cases{i}.rhs;
    expected = test_cases{i}.expected;
    actual = lhs.merge(rhs);
    assertLoop(expected == actual,i)
end

res = true;

% ------------------------------ END OF CODE ------------------------------
