function res = test_levelSet_display
% test_levelSet_display - unit test function of display
%
% Syntax:
%    res = test_levelSet_display
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

res = true;

% init symbolic variables
syms a b x y

% empty set
ls = levelSet.empty(2)

% single equation with different comparison operators
eq = x^2 + y^2 - 4;
ls = levelSet(eq,[x;y],'==')
ls = levelSet(eq,[x;y],'<')
ls = levelSet(eq,[x;y],'<=')

% multiple equations with different comparison operators
eq1 = sin(a) + log(b);
eq2 = abs(a)*b;
ls = levelSet([eq1;eq2],[a;b],{'<=';'<='})
ls = levelSet([eq1;eq2],[a;b],{'<=';'<'})

% independent variables
eq1 = a; eq2 = b; eq3 = x; eq4 = y;
ls = levelSet([eq1;eq2;eq3;eq4],[a;b;x;y],{'<';'<=';'<';'<='})

% unused variable in vars
eq = a + y;
ls = levelSet(eq,[a;b;y],'==')

% ------------------------------ END OF CODE ------------------------------
