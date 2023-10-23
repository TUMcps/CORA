function res = test_levelSet_ne
% test_levelSet_ne - unit test function of '~=' operator
%
% Syntax:
%    res = test_levelSet_ne
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
% Written:       28-April-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% assume true
res = true;

% init symbolic variables
syms a b c x y z

% init level set
eqs = -x^2 - y^2 + 5;
ls1 = levelSet(eqs,[x;y],'<=');

% different order in symbolic function
eqs = 5 - y^2 - x^2;
ls1_ = levelSet(eqs,[x;y],'<=');

% different variable names
eqs = -a^2 - b^2 + 5;
ls1__ = levelSet(eqs,[a;b],'<=');

% compare equal sets
res(end+1,1) = ~(ls1 ~= ls1_);
res(end+1,1) = ~ne(ls1,ls1__);


% different level set function
eqs = sin(a) + cos(b);
ls2 = levelSet(eqs,[a;b],'<=');

res(end+1,1) = ls1 ~= ls2;


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

res(end+1,1) = ~(ls1 ~= ls2);


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
