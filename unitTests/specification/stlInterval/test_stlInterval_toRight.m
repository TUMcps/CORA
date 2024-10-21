function res = test_stlInterval_toRight
% test_stlInterval_toRight - unit test function of toRight
%
% Syntax:
%    res = test_stlInterval_toRight
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
    {stlInterval(0,1,true,true), stlInterval(1,inf,false,false)};
    {stlInterval(0,1,true,false), stlInterval(1,inf,true,false)};
    {stlInterval(3,inf,true,false), stlInterval()};
};

% run tests
for i = 1:length(test_cases)
    int = test_cases{i}{1};
    expected = test_cases{i}{2};
    actual = toRight(int);
    assertLoop(actual == expected,i)
end

% empty interval
try
    toRight(stlInterval());
    assert(false);
catch ME
    if ~strcmp(ME.identifier,'CORA:emptySet')
        rethrow(ME)
    end
end

res = true;

% ------------------------------ END OF CODE ------------------------------
