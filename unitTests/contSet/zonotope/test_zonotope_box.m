function res = test_zonotope_box
% test_zonotope_box - unit test function of box
%
% Syntax:
%    res = test_zonotope_box
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
% Written:       26-August-2019
% Last update:   09-August-2020 (enhance randomness of test)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% set tolerance
tol = 1e-9;

% init zonotope
Z = zonotope([1;0],[2 -1; 4 1]);

% compute axis-aligned box
Zbox = box(Z);

% convert to interval and back to zonotope
Ztrue = zonotope([1;0],[3 0; 0 5]);

% check if axis-aligned box same as interval
res = compareMatrices([Zbox.c, Zbox.G],[Ztrue.c, Ztrue.G],tol);

% ------------------------------ END OF CODE ------------------------------
