function res = test_levelSet_isemptyobject
% test_levelSet_isemptyobject - unit test function of isemptyobject
%
% Syntax:
%    res = test_levelSet_isemptyobject
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
% See also: -

% Authors:       Mark Wetzlinger
% Written:       03-June-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true(0);

% empty level set
ls = levelSet.empty(2);
res(end+1,1) = ~isemptyobject(ls);

% 2D level set
syms x y
eq = x + y - 1;
ls2 = levelSet(eq,[x;y],'==');
res(end+1,1) = ~isemptyobject(ls2);

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
