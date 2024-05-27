function res = test_levelSet_and
% test_levelSet_and - unit test function of intersection
%
% Syntax:
%    res = test_levelSet_and
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
% See also: ------

% Authors:       Maximilian Perschl
% Written:       08-November-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;

% intersection of an interval with a linear levelSet
syms x y
eq = x + y - 1;
ls = levelSet(eq,[x;y],'==');
myInt = interval([0.5;0],[1.5;1]);
intSet = ls & myInt;
expectedSolution = interval([0.5;0],[1.0;0.5]);
if ~isequal(intSet,expectedSolution,1e-10)
    res = false;
end

% intersection of an interval with a non-linear levelSet
eq = x^2 + y - 1;
ls = levelSet(eq,[x;y],'==');
intSet = ls & myInt;
expectedSolution = interval([0.5;0],[1.0;0.8125]);
if ~isequal(intSet,expectedSolution,1e-10)
    res = false;
end

% intersection of two non-linear level sets
ls1 = levelSet([x^2 + y; y^2 - x],[x;y],'<=');
ls2 = levelSet(x^3 + y,[x;y],'<=');
lsRes = ls1 & ls2;
display(lsRes)                          % previously resulted in an error
if ~contains(lsRes,[0.5;-0.5]) || length(lsRes.eq) ~= length(lsRes.compOp)
    res = false;
end

% ------------------------------ END OF CODE ------------------------------
