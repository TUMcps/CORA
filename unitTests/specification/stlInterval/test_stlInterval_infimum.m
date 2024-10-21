function res = test_stlInterval_infimum
% test_stlInterval_infimum - unit test function of infimum
%
% Syntax:
%    res = test_stlInterval_infimum
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
    % {int, infi, isMin}
    {stlInterval(0,1,true,true), 0, true};
    {stlInterval(0,1,false,true), 0, false};
    {stlInterval(0,1,true,false), 0, true};
    {stlInterval(0,1,false,false), 0, false};
};

% run tests
for i = 1:length(test_cases)
    int = test_cases{i}{1};
    expectedInfi = test_cases{i}{2};
    expectedIsMin = test_cases{i}{3};
    [infi,isMin] = infimum(int);
    assertLoop(withinTol(infi,expectedInfi,eps),i)
    assertLoop(isMin == expectedIsMin,i)
end

% empty interval
assert(isempty(infimum(stlInterval())))

res = true;

% ------------------------------ END OF CODE ------------------------------
