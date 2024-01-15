function res = test_conZonotope_empty
% test_conZonotope_empty - unit test function of empty
%
% Syntax:
%    res = test_conZonotope_empty
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
% Written:       09-January-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true(0);

% empty case
cZ = conZonotope.empty(3);
res(end+1,1) = dim(cZ) == 3;
res(end+1,1) = representsa(cZ,'emptySet');

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
