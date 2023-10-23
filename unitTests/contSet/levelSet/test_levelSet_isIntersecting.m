function res = test_levelSet_isIntersecting
% test_levelSet_isIntersecting - unit test function of isIntersecting
%
% Syntax:
%    res = test_levelSet_isIntersecting
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
% See also: isIntersecting

% Authors:       Maximilian Perschl
% Written:       08-November-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;
% Define problem
% intersection with linear inequality levelSet
syms x y
eq = x + y - 1;
ls = levelSet(eq,[x;y],'<=');
myInt = interval([0.5;0],[1.5;1]);
if ~isIntersecting(ls,myInt,"approx")
    res = false;
end
% intersection with linear equality levelSet
ls = levelSet(eq,[x;y],'==');
if ~isIntersecting(ls,myInt,"approx")
    res = false;
end
% different interval, not intersecting:
myInt = interval([1+1e-5;0],[1.2;1]);
if isIntersecting(ls,myInt,"approx")
    res = false;
end

% tests for nonlinear inequality levelSet
eq = 3*x + y^2 - 1;
ls = levelSet(eq,[x;y],'<=');
myInt = interval([0.25;0.5],[0.5;1]);
if ~isIntersecting(ls,myInt,"approx")
    res = false;
end
% different interval, not intersecting:
myInt = interval([0.25+1e-3;0.5],[1;1]);
if isIntersecting(ls,myInt,"approx")
    res = false;
end


% ------------------------------ END OF CODE ------------------------------
