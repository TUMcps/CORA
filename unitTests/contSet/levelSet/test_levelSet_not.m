function res = test_levelSet_not
% test_levelSet_not - unit test function of complement operation
%
% Syntax:
%    res = test_levelSet_not
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
% Written:       25-April-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% assume true
res = true;

% init symbolic variables
syms a b x y

% init level set (<=)
eq = -x^2 - y^2 + 5;
ls = levelSet(eq,[x;y],'<=');

% compute complement
ls_ = ~ls;

% true solution
eq = x^2 + y^2 - 5;
ls_true = levelSet(eq,[x;y],'<');

% check
res(end+1,1) = isequal(ls_,ls_true);


% another level set (<)
eq = sin(a) + cos(b);
ls = levelSet(eq,[a;b],'<');

% compute complement
ls_ = ~ls;

% true solution
eq = -sin(a) - cos(b);
ls_true = levelSet(eq,[a;b],'<=');

% check
res(end+1,1) = isequal(ls_,ls_true);


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
