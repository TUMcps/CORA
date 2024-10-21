function res = test_stlInterval_and
% test_stlInterval_and - unit test function of intersection
%
% Syntax:
%    res = test_stlInterval_and
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
    % {a, b, expected}
    {stlInterval(0,1,true,true), stlInterval(0,1,true,true), stlInterval(0,1,true,true)};
    {stlInterval(0,1,true,true), stlInterval(0,0.5,true,true), stlInterval(0,0.5,true,true)};
    {stlInterval(0,1,true,true), stlInterval(0.5,1,true,true), stlInterval(0.5,1,true,true)};
    {stlInterval(0,1,true,true), stlInterval(0,1,false,true), stlInterval(0,1,false,true)};
    {stlInterval(0,1,true,true), stlInterval(0,1,true,false), stlInterval(0,1,true,false)};
    {stlInterval(0,1,true,true), stlInterval(0,0.5,false,true), stlInterval(0,0.5,false,true)};
    {stlInterval(0,1,true,true), stlInterval(0.5,1,true,false), stlInterval(0.5,1,true,false)};
    {stlInterval(0,1,true,true), stlInterval(0,1,false,false), stlInterval(0,1,false,false)};
    {stlInterval(0,1,true,true), stlInterval(1,2,true,true), stlInterval(1,1,true,true)};
    {stlInterval(0,1,true,true), stlInterval(1,2,false,true), stlInterval()};
    {stlInterval(0,1,true,true), stlInterval(), stlInterval()};
    {stlInterval(0,1,true,true), stlInterval(0,inf,true,false), stlInterval(0,1,true,true)};
};

% run tests
for i = 1:length(test_cases)
    a = test_cases{i}{1};
    b = test_cases{i}{2};
    expected = test_cases{i}{3};
    actual = a & b;
    assertLoop(actual == expected,i)
end

res = true;

% ------------------------------ END OF CODE ------------------------------
