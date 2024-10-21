function res = test_kleeneSignal_findIntervals
% test_kleeneSignal_findIntervals - unit test function of findIntervals
%
% Syntax:
%    res = test_kleeneSignal_findIntervals
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
tt = true;
ff = false;
boolSig = pointSegmentSignal([0,1,2,2.5,3,4],[tt,tt,ff,tt,tt,ff,ff,tt,ff,ff,ff,tt]);
isUnknown = pointSegmentSignal.indicator(stlInterval(1.5,2.25,true,false),true,false);
sig = kleeneSignal(boolSig,isUnknown);

% test case definition
test_cases = {
    % {sig, match, expected}
    {sig, kleene.True, [stlInterval(0,1,true,false),stlInterval(1,1.5,false,false),stlInterval(2.5,3,false),stlInterval(4,inf,false)]};
    {sig, kleene.Unknown, [stlInterval(1.5,2.25,true,false)]};
    {sig, kleene.False, [stlInterval(1),stlInterval(2.25,2.5,true,true),stlInterval(3,4,true)]};
};

% run tests
for i = 1:length(test_cases)
    sig = test_cases{i}{1};
    match = test_cases{i}{2};
    expected = test_cases{i}{3};
    actual = sig.findIntervals(match);
    assertLoop(compareIntervals(actual,expected),i)
end

res = true;

% ------------------------------ END OF CODE ------------------------------
