function res = test_kleeneSignal_at
% test_kleeneSignal_at - unit test function of at
%
% Syntax:
%    res = test_kleeneSignal_at
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
boolSig = pointSegmentSignal([0,1,2,2.5],[true,true,false,true,true,false,false,true]);
isUnknown = pointSegmentSignal.indicator(stlInterval(1.5,2.25,true,false),true,false);
sig = kleeneSignal(boolSig,isUnknown);

% test case definition
tt = kleene.True;
uu = kleene.Unknown;
ff = kleene.False;
test_cases = {
    % {time, expected}
    {0, tt};
    {0.5, tt};
    {1, ff};
    {1.25, tt};
    {1.5, uu};
    {2, uu};
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
