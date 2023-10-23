function res = test_levelSet_isequal
% test_levelSet_isequal - unit test function of equality check
%
% Syntax:
%    res = test_levelSet_isequal
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

% Authors:       Mark Wetzlinger
% Written:       11-January-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% assume true
res = true;

% init symbolic variables
syms a b c x y z

% init level set
eq = -x^2 - y^2 + 5;
ls1 = levelSet(eq,[x;y],'<=');

% different order in symbolic function
eq = 5 - y^2 - x^2;
ls1_ = levelSet(eq,[x;y],'<=');

% different variable names
eq = -a^2 - b^2 + 5;
ls1__ = levelSet(eq,[a;b],'<=');

% compare equal sets
if ~isequal(ls1,ls1_)
    res = false;
elseif ~isequal(ls1,ls1__)
    res = false;
end


% different level set function
eq = sin(a) + cos(b);
ls2 = levelSet(eq,[a;b],'<=');

if isequal(ls1,ls2)
    res = false;
end


% multiple equations
eq1 = -x^2 + y^2 - 5;
eq2 = sin(x) + cos(y);
eq3 = 2*z;
ls1 = levelSet([eq1;eq2;eq3],[x;y;z],{'<=','<','<='});

% different order and variable names
eq1 = cos(b) + sin(a);
eq2 = 2*c;
eq3 = b^2 - 5 - a^2;
ls2 = levelSet([eq1;eq2;eq3],[a;b;c],{'<','<=','<='});

if ~isequal(ls1,ls2)
    res = false;
end


% ------------------------------ END OF CODE ------------------------------
