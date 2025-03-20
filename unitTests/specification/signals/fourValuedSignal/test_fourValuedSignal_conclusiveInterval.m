function res = test_fourValuedSignal_conclusiveInterval
% test_fourValuedSignal_conclusiveInterval - unit test function of conclusiveInterval
%
% Syntax:
%    res = test_fourValuedSignal_conclusiveInterval
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

test_cases = {
    % {sig, ignore, expected}
    {fourValuedSignal.indicator(stlInterval(0,1),tt,ii), stlInterval(), stlInterval(0,1)};
    {fourValuedSignal.indicator(stlInterval(0,1),uu,ii), stlInterval(), stlInterval(0,1)};
    {fourValuedSignal.indicator(stlInterval(0,1),ff,ii), stlInterval(), stlInterval(0,1)};
    {fourValuedSignal.indicator(stlInterval(0,1),ii,ii), stlInterval(), stlInterval()};
    {fourValuedSignal.indicator(stlInterval(0,1,false),tt,ii), stlInterval(), stlInterval()};
    {fourValuedSignal.indicator(stlInterval(0,1,false),tt,ii), stlInterval(0), stlInterval(0,1,true,false)};
    {fourValuedSignal.indicator(stlInterval(0,1),tt,ii).set(stlInterval(4,5),ff), stlInterval(), stlInterval(0,1)};
    {fourValuedSignal.indicator(stlInterval(0,1),tt,ii).set(stlInterval(4,5),ff), stlInterval(1,3), stlInterval(0,3)};
    {fourValuedSignal.indicator(stlInterval(0,1),tt,ii).set(stlInterval(4,5),ff), stlInterval(2,3), stlInterval(0,1)};
    {fourValuedSignal.indicator(stlInterval(0,1),tt,ii).set(stlInterval(4,5),ff), stlInterval(1,4,false), stlInterval(0,5)};
    {fourValuedSignal.indicator(stlInterval(0,1,true,false),tt,ii).set(stlInterval(4,5,false),ff), stlInterval(1,4,false), stlInterval(0,1,true,false)};
};

% run test cases
for i = 1:length(test_cases)
    sig = test_cases{i}{1};
    ignore = test_cases{i}{2};
    expected = test_cases{i}{3};
    actual = sig.conclusiveInterval(ignore);
    assertLoop(expected == actual,i)
end

res = true;

% ------------------------------ END OF CODE ------------------------------
