function res = test_halfspace_isemptyobject
% test_halfspace_isemptyobject - unit test function of isemptyobject
%
% Syntax:
%    res = test_halfspace_isemptyobject
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

% empty halfspace
hs = halfspace.empty(2);
res(end+1,1) = ~isemptyobject(hs);

% 3D halfspace
a = [3; 2; -1];
b = 0.5;
hs = halfspace(a,b);
res(end+1,1) = ~isemptyobject(hs);

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
