function res = test_stlInterval_rightClosure
% test_stlInterval_rightClosure - unit test function of rightClosure
%
% Syntax:
%    res = test_stlInterval_rightClosure
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
    {stlInterval(0,1,false,true), stlInterval(0,1,false,true)};
    {stlInterval(0,1,true,false), stlInterval(0,1,true,true)};
    {stlInterval(0,1,false,false), stlInterval(0,1,false,true)};
    {stlInterval(), stlInterval()};
    {stlInterval(0,inf), stlInterval(0,inf)};
};

% run tests
for i = 1:length(test_cases)
    int = test_cases{i}{1};
    expected = test_cases{i}{2};
    actual = rightClosure(int);
    assertLoop(actual == expected,i)
end

res = true;

% ------------------------------ END OF CODE ------------------------------
