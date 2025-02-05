function res = test_stlInterval_contains
% test_stlInterval_contains - unit test function of contains
%
% Syntax:
%    res = test_stlInterval_contains
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

tol = 1e-8;

% test case definition
test_cases_exact = {
    % {a, b, expected}
    % combinations of open/closed intervals
    {stlInterval(0,1,true,true), stlInterval(0,1,true,true), true};
    {stlInterval(0,1,true,true), stlInterval(0,1,true,false), true};
    {stlInterval(0,1,true,true), stlInterval(0,1,false,true), true};
    {stlInterval(0,1,true,true), stlInterval(0,1,false,false), true};
    {stlInterval(0,1,true,false), stlInterval(0,1,true,true), false};
    {stlInterval(0,1,false,true), stlInterval(0,1,true,true), false};
    {stlInterval(0,1,true,true), stlInterval(0,0.5,true,true), true};
    
    % containment of CORA intervals
    {stlInterval(0,1,true,true), interval(0,1), true};
    {stlInterval(0,1,true,true), interval(0.5,0.75), true};

    % containment of points
    {stlInterval(0,1,true,true), 0.5, true};
    {stlInterval(0,1,true,true), 1.5, false};
    {stlInterval(0,1,true,true), 0, true};
    {stlInterval(0,1,true,true), 1, true};
    {stlInterval(0,1,false,false), 0.5, true};
    {stlInterval(0,1,false,false), 1.5, false};
    {stlInterval(0,1,false,false), 0, false};
    {stlInterval(0,1,false,false), 1, false};
};

% define test cases for approximate containment
test_cases_approx = {
    {stlInterval(0,1,false,false), interval(0,1), false};
    {stlInterval(0,1,false,false), interval(0.5,0.75), true};
};

% run exact test cases
for i = 1:length(test_cases_exact)
    a = test_cases_exact{i}{1};
    b = test_cases_exact{i}{2};
    expected = test_cases_exact{i}{3};
    actual = contains(a,b,'exact',tol);
    assertLoop(actual == expected,i);
end

% run approximate test cases
for i = 1:length(test_cases_approx)
    a = test_cases_approx{i}{1};
    b = test_cases_approx{i}{2};
    expected = test_cases_approx{i}{3};
    actual = contains(a,b,'approx',tol);
    assertLoop(actual == expected,i);
end

res = true;

% ------------------------------ END OF CODE ------------------------------
