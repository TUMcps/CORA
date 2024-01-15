function res = test_zonoBundle_volume
% test_zonoBundle_volume - unit test function of volume
%
% Syntax:
%    res = test_zonoBundle_volume
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
% Written:       23-April-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% fully-empty zonoBundle
zB = zonoBundle.empty(2);
val = volume(zB);
res = val == 0;

% non-empty intersection
Z1 = zonotope([1;1], [3 0; 0 2]);
Z2 = zonotope([0;0], [2 2; 2 -2]);
zB = zonoBundle({Z1,Z2});
val = volume(zB);
% true result computed by hand
val_true = 6*4 - 1*1/2 - 1*1/2 - 3*3/2;
res(end+1,1) = withinTol(val,val_true);

% empty intersection
Z2 = zonotope([-4;1],[0.5 1; 1 -1]);
zB = zonoBundle({Z1,Z2});
val = volume(zB);
res(end+1,1) = val == 0;

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
