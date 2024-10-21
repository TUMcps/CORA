function res = test_fourValuedSignal_at
% test_fourValuedSignal_at - unit test function of at
%
% Syntax:
%    res = test_fourValuedSignal_at
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
boolSig = pointSegmentSignal([0,1,2,2.5],[true,true,false,true,true,false,false,true]);
isUnknown = pointSegmentSignal.indicator(stlInterval(1.5,2.25,true,false),true,false);
isInconclusive = pointSegmentSignal.indicator(stlInterval(1.75,2.1,false,true),true,false);
sig = fourValuedSignal(kleeneSignal(boolSig,isUnknown),isInconclusive);

% test case definition
tt = fourValued.True;
uu = fourValued.Unknown;
ff = fourValued.False;
ii = fourValued.Inconclusive;
test_cases = {
    % {time, expected}
    {0, tt};
    {0.5, tt};
    {1, ff};
    {1.25, tt};
    {1.5, uu};
    {1.75, uu};
    {2, ii};
    {2.1, ii};
    {2.25, ff};
    {2.5, ff};
    {3, tt};
    {inf, tt};
};

% run tests
for i = 1:length(test_cases)
    time = test_cases{i}{1};
    expected = test_cases{i}{2};
    actual = sig.at(time);
    assertLoop(actual == expected,i)
end

res = true;

% ------------------------------ END OF CODE ------------------------------
