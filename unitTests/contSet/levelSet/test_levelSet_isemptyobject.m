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

% Author:       Mark Wetzlinger
% Written:      03-June-2022
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% instantiate level sets
ls1 = levelSet();

syms x y
eq = x + y - 1;
ls2 = levelSet(eq,[x;y],'==');

% check results
res = isemptyobject(ls1) && ~isemptyobject(ls2);

%------------- END OF CODE --------------
