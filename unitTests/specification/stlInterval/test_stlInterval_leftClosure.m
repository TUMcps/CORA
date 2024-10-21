function res = test_stlInterval_leftClosure
% test_stlInterval_leftClosure - unit test function of leftClosure
%
% Syntax:
%    res = test_stlInterval_leftClosure
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
    {stlInterval(0,1,true,true), stlInterval(0,1,true,true)};
    {stlInterval(0,1,false,true), stlInterval(0,1,true,true)};
    {stlInterval(0,1,true,false), stlInterval(0,1,true,false)};
    {stlInterval(0,1,false,false), stlInterval(0,1,true,false)};
    {stlInterval(), stlInterval()};
};

% run tests
for i = 1:length(test_cases)
    int = test_cases{i}{1};
    expected = test_cases{i}{2};
    actual = leftClosure(int);
    assertLoop(actual == expected,i)
end

res = true;

% ------------------------------ END OF CODE ------------------------------
