function res = test_levelSet_copy
% test_levelSet_copy - unit test function of copy
%
% Syntax:
%    res = test_levelSet_copy
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
% Written:       02-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% 2D levelSet
syms x y
eq = x^2 + y^2 - 4;
ls = levelSet(eq,[x;y],'==');
ls_copy = copy(ls);
assert(isequal(ls,ls_copy));


% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
