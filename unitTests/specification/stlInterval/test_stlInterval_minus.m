function res = test_stlInterval_minus
% test_stlInterval_minus - unit test function of minus
%
% Syntax:
%    res = test_stlInterval_minus
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
    {stlInterval(10,20,true,true), stlInterval(3,4,true,true), stlInterval(6,17,true,true)};
    {stlInterval(10,20,true,true), stlInterval(3,4,false,true), stlInterval(6,17,true,false)};
    {stlInterval(10,20,true,true), stlInterval(3,4,true,false), stlInterval(6,17,false,true)};
    {stlInterval(10,20,true,true), stlInterval(3,4,false,false), stlInterval(6,17,false,false)};

    {stlInterval(10,20,true,false), stlInterval(3,4,true,true), stlInterval(6,17,true,false)};
    {stlInterval(10,20,true,false), stlInterval(3,4,false,true), stlInterval(6,17,true,false)};
    {stlInterval(10,20,true,false), stlInterval(3,4,true,false), stlInterval(6,17,false,false)};
    {stlInterval(10,20,true,false), stlInterval(3,4,false,false), stlInterval(6,17,false,false)};

    {stlInterval(10,20,false,true), stlInterval(3,4,true,true), stlInterval(6,17,false,true)};
    {stlInterval(10,20,false,true), stlInterval(3,4,false,true), stlInterval(6,17,false,false)};
    {stlInterval(10,20,false,true), stlInterval(3,4,true,false), stlInterval(6,17,false,true)};
    {stlInterval(10,20,false,true), stlInterval(3,4,false,false), stlInterval(6,17,false,false)};

    {stlInterval(10,20,false,false), stlInterval(3,4,true,true), stlInterval(6,17,false,false)};
    {stlInterval(10,20,false,false), stlInterval(3,4,false,true), stlInterval(6,17,false,false)};
    {stlInterval(10,20,false,false), stlInterval(3,4,true,false), stlInterval(6,17,false,false)};
    {stlInterval(10,20,false,false), stlInterval(3,4,false,false), stlInterval(6,17,false,false)};

    {stlInterval(10,20), stlInterval(25,30), stlInterval()};
    {stlInterval(10,20), stlInterval(20,30), stlInterval(0)};
    {stlInterval(10,20), stlInterval(15,30), stlInterval(0,5)};
    {stlInterval(10,20), stlInterval(3,inf), stlInterval(0,17)};

    {stlInterval(), stlInterval(3,4), stlInterval()};
    {stlInterval(3,4), stlInterval(), stlInterval()};

    {stlInterval(10,20), 3, stlInterval(7,17)};
    {stlInterval(10,20), interval(3,4), stlInterval(6,17)};
};

% run tests
for i = 1:length(test_cases)
    a = test_cases{i}{1};
    b = test_cases{i}{2};
    expected = test_cases{i}{3};
    actual = a - b;
    assertLoop(actual == expected,i)
end

% test exception
try
    stlInterval(10,20) - interval([3 4], [5 6]); %#ok<VUNUS>
    assert(false);
catch ME
    if ~strcmp(ME.identifier, 'CORA:dimensionMismatch')
        rethrow(ME)
    end
end

res = true;

% ------------------------------ END OF CODE ------------------------------
