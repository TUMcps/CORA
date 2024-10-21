function res = test_stlInterval_toLeft
% test_stlInterval_toLeft - unit test function of toLeft
%
% Syntax:
%    res = test_stlInterval_toLeft
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
    {stlInterval(0,1,true,true), stlInterval()};
    {stlInterval(0,1,false,true), stlInterval(0)};
    {stlInterval(3,5,true,false), stlInterval(0,3,true,false)};
    {stlInterval(3,5,false,false), stlInterval(0,3,true,true)};
};

% run tests
for i = 1:length(test_cases)
    int = test_cases{i}{1};
    expected = test_cases{i}{2};
    actual = toLeft(int);
    assertLoop(actual == expected,i)
end

% empty interval
try
    toLeft(stlInterval());
    assert(false);
catch ME
    if ~strcmp(ME.identifier,'CORA:emptySet')
        rethrow(ME)
    end
end

res = true;

% ------------------------------ END OF CODE ------------------------------
