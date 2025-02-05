function res = test_stlInterval_minkDiff
% test_stlInterval_minkDiff - unit test function of minkDiff
%
% Syntax:
%    res = test_stlInterval_minkDiff
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
    
    % lhs closed/closed vs all possibilities for rhs
    {stlInterval(10,20,true,true), stlInterval(3,4,true,true), stlInterval(7,16,true,true)};
    {stlInterval(10,20,true,true), stlInterval(3,4,false,true), stlInterval(7,16,true,true)};
    {stlInterval(10,20,true,true), stlInterval(3,4,true,false), stlInterval(7,16,true,true)};
    {stlInterval(10,20,true,true), stlInterval(3,4,false,false), stlInterval(7,16,true,true)};

    % lhs closed/open vs all possibilities for rhs
    {stlInterval(10,20,true,false), stlInterval(3,4,true,true), stlInterval(7,16,true,false)};
    {stlInterval(10,20,true,false), stlInterval(3,4,false,true), stlInterval(7,16,true,false)};
    {stlInterval(10,20,true,false), stlInterval(3,4,true,false), stlInterval(7,16,true,true)};
    {stlInterval(10,20,true,false), stlInterval(3,4,false,false), stlInterval(7,16,true,true)};

    % lhs open/closed vs all possibilities for rhs
    {stlInterval(10,20,false,true), stlInterval(3,4,true,true), stlInterval(7,16,false,true)};
    {stlInterval(10,20,false,true), stlInterval(3,4,false,true), stlInterval(7,16,true,true)};
    {stlInterval(10,20,false,true), stlInterval(3,4,true,false), stlInterval(7,16,false,true)};
    {stlInterval(10,20,false,true), stlInterval(3,4,false,false), stlInterval(7,16,true,true)};

    % lhs open/open vs all possibilities for rhs
    {stlInterval(10,20,false,false), stlInterval(3,4,true,true), stlInterval(7,16,false,false)};
    {stlInterval(10,20,false,false), stlInterval(3,4,false,true), stlInterval(7,16,true,false)};
    {stlInterval(10,20,false,false), stlInterval(3,4,true,false), stlInterval(7,16,false,true)};
    {stlInterval(10,20,false,false), stlInterval(3,4,false,false), stlInterval(7,16,true,true)};

    % empty result
    {stlInterval(10,20), stlInterval(20,30), stlInterval()};
    {stlInterval(10,20), stlInterval(15,30), stlInterval()};
    {stlInterval(10,20), stlInterval(3,inf), stlInterval()};
    {stlInterval(), stlInterval(3,4), stlInterval()};

    {stlInterval(10,20), stlInterval(15,17), stlInterval(0,3)};
    {stlInterval(3,4), stlInterval(), stlInterval(0,inf)};
};

% run tests
for i = 1:length(test_cases)
    a = test_cases{i}{1};
    b = test_cases{i}{2};
    expected = test_cases{i}{3};
    actual = minkDiff(a,b);
    assertLoop(actual == expected,i)
end

res = true;

% ------------------------------ END OF CODE ------------------------------
