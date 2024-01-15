function res = test_levelSet_lift
% test_levelSet_lift - unit test function of lift
%
% Syntax:
%    res = test_levelSet_lift
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

% init symbolic variables
syms a b c d
syms s1 s2 s4 s5 s6 s8

% init 2D level set
eq = a^2 + b^2 - 4;
ls = levelSet(eq,[a;b],'<=');

% project to higher dimension
N = 4;
dims = [2 4];
ls_ = lift(ls,N,dims);

% true solution
eq_true = b^2 + d^2 - 4;
ls_true = levelSet(eq_true,[a;b;c;d],'<=');

% compare
res = isequal(ls_,ls_true);


% multiple equations
eq1 = a^2 - b;
eq2 = c + sqrt(d);
ls = levelSet([eq1;eq2],[a;b;c;d],{'<=','<'});

% project to higher dimension
N = 8;
dims = [2 3 5 7];
ls_ = lift(ls,N,dims);

% true solution
ls_true = levelSet([eq1;eq2],[s1;a;b;s4;c;s6;d;s8],{'<=','<'});

% compare
res(end+1,1) = isequal(ls_,ls_true);


% unused variable
eq = a^2 - c;
ls = levelSet(eq,[a;b;c],'==');

% project to higher dimension
N = 5;
dims = [1 3 4];
ls_ = lift(ls,N,dims);

% true solution
ls_true = levelSet(eq,[a;s2;b;c;s5],'==');

% compare
res(end+1,1) = isequal(ls_,ls_true);


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
