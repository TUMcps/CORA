function res = test_stlInterval_interval
% test_stlInterval_interval - unit test function of interval
%
% Syntax:
%    res = test_stlInterval_interval
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
    {stlInterval(0,1,true,true), interval(0,1)};
    {stlInterval(0,1,true,false), interval(0,1)};
    {stlInterval(0,1,false,true), interval(0,1)};
    {stlInterval(0,1,false,false), interval(0,1)};
    {stlInterval(0,inf), interval(0,inf)};
};

% run tests
for i = 1:length(test_cases)
    int = test_cases{i}{1};
    expected = test_cases{i}{2};
    actual = interval(int);
    assertLoop(actual == expected,i)
end

res = true;

% ------------------------------ END OF CODE ------------------------------
