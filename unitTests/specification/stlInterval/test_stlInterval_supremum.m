function res = test_stlInterval_supremum
% test_stlInterval_supremum - unit test function of supremum
%
% Syntax:
%    res = test_stlInterval_supremum
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

% test case definition
test_cases = {
    % {int, sup, isMax}
    {stlInterval(0,1,true,true), 1, true};
    {stlInterval(0,1,false,true), 1, true};
    {stlInterval(0,1,true,false), 1, false};
    {stlInterval(0,1,false,false), 1, false};
    {stlInterval(0,inf), inf, false};
};

% run tests
for i = 1:length(test_cases)
    int = test_cases{i}{1};
    expectedSup = test_cases{i}{2};
    expectedIsMax = test_cases{i}{3};
    [sup,isMax] = supremum(int);
    assertLoop(withinTol(sup,expectedSup,eps),i)
    assertLoop(isMax == expectedIsMax,i)
end

% empty interval
assert(isempty(supremum(stlInterval())))

res = true;

% ------------------------------ END OF CODE ------------------------------
