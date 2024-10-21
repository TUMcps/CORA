function res = test_fourValuedSignal_indicator
% test_fourValuedSignal_indicator - unit test function of indicator
%
% Syntax:
%    res = test_fourValuedSignal_indicator
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

% test case definition
tt = kleene.True;
uu = kleene.Unknown;
ff = kleene.False;
interval = stlInterval(1,2);
blank = pointSegmentSignal.indicator(stlInterval(),false,false);
inside = pointSegmentSignal.indicator(interval,true,false);
test_cases = {
    % {val, default, expected}
    {tt, tt, kleeneSignal(~blank,blank)};
    {tt, uu, kleeneSignal(inside,~inside)};
    {tt, ff, kleeneSignal(inside,blank)};
    {uu, tt, kleeneSignal(~inside,inside)};
    {uu, uu, kleeneSignal(blank,~blank)};
    {uu, ff, kleeneSignal(blank,inside)};
    {ff, tt, kleeneSignal(~inside,blank)};
    {ff, uu, kleeneSignal(blank,~inside)};
    {ff, ff, kleeneSignal(blank,blank)};
};

% run tests
for i = 1:length(test_cases)
    val = test_cases{i}{1};
    default = test_cases{i}{2};
    expected = test_cases{i}{3};
    actual = kleeneSignal.indicator(interval,val,default);
    assertLoop(actual == expected,i)
end

res = true;

% ------------------------------ END OF CODE ------------------------------
