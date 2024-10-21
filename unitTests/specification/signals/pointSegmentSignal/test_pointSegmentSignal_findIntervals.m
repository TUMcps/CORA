function res = test_pointSegmentSignal_findIntervals
% test_pointSegmentSignal_findIntervals - unit test function of findIntervals
%
% Syntax:
%    res = test_pointSegmentSignal_findIntervals
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
% Written:       16-February-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% signals
tt = true;
ff = false;
sig1 = pointSegmentSignal([0,1,2,2.5,3,4],[tt,tt,ff,tt,tt,ff,ff,tt,ff,ff,ff,tt]);
sig2 = ~sig1;

% test case definition
test_cases = {
    % {sig, cond, expected}
    {sig1, @(x) x, [stlInterval(0,1,true,false),stlInterval(1,2,false,true),stlInterval(2.5,3,false),stlInterval(4,inf,false)]};
    {sig1, @(x) ~x, [stlInterval(1),stlInterval(2,2.5,false,true),stlInterval(3,4,true)]};
    {sig2, @(x) x, [stlInterval(1),stlInterval(2,2.5,false,true),stlInterval(3,4,true)]};
    {sig2, @(x) ~x, [stlInterval(0,1,true,false),stlInterval(1,2,false,true),stlInterval(2.5,3,false),stlInterval(4,inf,false)]};
};

% run tests
for i = 1:length(test_cases)
    sig = test_cases{i}{1};
    cond = test_cases{i}{2};
    expected = test_cases{i}{3};
    actual = sig.findIntervals(cond);
    assertLoop(all(arrayfun(@(a,e) a == e,actual,expected)),i)
end

res = true;

% ------------------------------ END OF CODE ------------------------------
