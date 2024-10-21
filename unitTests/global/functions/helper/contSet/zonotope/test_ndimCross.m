function res = test_ndimCross()
% test_ndimCross - unit test function for n-dimensional cross product
%
% Syntax:
%    res = test_ndimCross()
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
% Written:       01-August-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;

% 2D
M = [5; -2];
v = ndimCross(M);
v_true = [-2; -5];
assert(all(withinTol(v,v_true)));

% 3D
M = [1 2; 3 4; 5 6];
v = ndimCross(M);
v_true = [-2; 4; -2];
assert(all(withinTol(v,v_true)));


% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
