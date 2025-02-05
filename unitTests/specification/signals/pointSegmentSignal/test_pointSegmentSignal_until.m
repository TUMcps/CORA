function res = test_pointSegmentSignal_until
% test_pointSegmentSignal_until - unit test function of until
%
% Syntax:
%    res = test_pointSegmentSignal_until
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
test_cases = {};

% case 1
test_cases{end+1}.lhs = pointSegmentSignal([0,1,3],[tt,tt,ff,ff,tt,tt]);
test_cases{end}.int = stlInterval(0,1,true);
test_cases{end}.rhs = pointSegmentSignal([0,4,6],[tt,tt,ff,ff,tt,tt]);
test_cases{end}.exp = pointSegmentSignal([0,4,5],[tt,tt,ff,ff,tt,tt]);

% case 2
test_cases{end+1}.lhs = pointSegmentSignal([0,1,3],[tt,tt,ff,ff,tt,tt]);
test_cases{end}.int = stlInterval(0,1,false);
test_cases{end}.rhs = pointSegmentSignal([0,4,6],[tt,tt,ff,ff,tt,tt]);
test_cases{end}.exp = pointSegmentSignal([0,1,3,4,5],[tt,tt,ff,ff,tt,tt,ff,ff,ff,tt]);

% case 3
test_cases{end+1}.lhs = pointSegmentSignal([0,2],[ff,ff,tt,tt]);
test_cases{end}.int = stlInterval(0,1,true);
test_cases{end}.rhs = pointSegmentSignal([0,2],[tt,tt,ff,ff]);
test_cases{end}.exp = pointSegmentSignal([0,2],[tt,tt,ff,ff]);

% case 4
test_cases{end+1}.lhs = pointSegmentSignal([0,2],[ff,ff,tt,tt]);
test_cases{end}.int = stlInterval(0,1,false);
test_cases{end}.rhs = pointSegmentSignal([0,2],[tt,tt,ff,ff]);
test_cases{end}.exp = pointSegmentSignal(0,[ff,ff]);

% case 5
test_cases{end+1}.lhs = pointSegmentSignal([0,2],[ff,tt,ff,ff]);
test_cases{end}.int = stlInterval(0,3,true);
test_cases{end}.rhs = pointSegmentSignal([0,2],[ff,ff,tt,tt]);
test_cases{end}.exp = pointSegmentSignal(0,[tt,tt]);

% case 6
test_cases{end+1}.lhs = pointSegmentSignal([0,2],[ff,ff,tt,tt]);
test_cases{end}.int = stlInterval(1,3,true);
test_cases{end}.rhs = pointSegmentSignal([0,2,5],[tt,tt,ff,ff,ff,tt]);
test_cases{end}.exp = pointSegmentSignal([0,2],[ff,ff,ff,tt]);

% case 7
test_cases{end+1}.lhs = pointSegmentSignal([0,2],[ff,ff,tt,tt]);
test_cases{end}.int = stlInterval(1,3,true);
test_cases{end}.rhs = pointSegmentSignal([0,2,5],[tt,tt,ff,ff,tt,tt]);
test_cases{end}.exp = pointSegmentSignal([0,2],[ff,ff,tt,tt]);

% case 8
test_cases{end+1}.lhs = pointSegmentSignal([0,2],[ff,ff,tt,tt]);
test_cases{end}.int = stlInterval(1,3,true,false);
test_cases{end}.rhs = pointSegmentSignal([0,2,5],[tt,tt,ff,ff,tt,tt]);
test_cases{end}.exp = pointSegmentSignal([0,2],[ff,ff,ff,tt]);

% run tests
for i = 1:length(test_cases)
    c = test_cases{i};
    expected = c.exp;
    actual = until(c.lhs,c.int,c.rhs);
    assertLoop(actual == expected,i)
end

res = true;

% ------------------------------ END OF CODE ------------------------------
