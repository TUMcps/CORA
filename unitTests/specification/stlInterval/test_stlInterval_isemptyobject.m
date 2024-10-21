function res = test_stlInterval_isemptyobject
% test_stlInterval_isemptyobject - unit test function of isemptyobject
%
% Syntax:
%    res = test_stlInterval_isemptyobject
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
    % {int, expected}
    {stlInterval(0,1,true,true), false};
    {stlInterval(0,1,true,false), false};
    {stlInterval(0,1,false,true), false};
    {stlInterval(0,1,false,false), false};
    {stlInterval(), true};
    {stlInterval(0,0,false,false), true};
};

% run tests
for i = 1:length(test_cases)
    int = test_cases{i}{1};
    expected = test_cases{i}{2};
    actual = isemptyobject(int);
    assertLoop(actual == expected,i)
end

res = true;

% ------------------------------ END OF CODE ------------------------------
