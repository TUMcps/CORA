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

% empty case
cZ = conZonotope.empty(3);
assert(dim(cZ) == 3);
assert(representsa(cZ,'emptySet'));

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
