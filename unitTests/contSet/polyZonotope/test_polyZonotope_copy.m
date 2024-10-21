function res = test_polyZonotope_copy
% test_polyZonotope_copy - unit test function of copy
%
% Syntax:
%    res = test_polyZonotope_copy
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

% 2D polyZonotope
c = [-1;3];
G = [-1 -1 3 2; 2 -1 -1 0];
E = [1 0 1 5; 0 1 3 0];
GI = [2; -1];
pZ = polyZonotope(c,G,GI,E);
pZ_copy = copy(pZ);
assert(isequal(pZ,pZ_copy));


% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
