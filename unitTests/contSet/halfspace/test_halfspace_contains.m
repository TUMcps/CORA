function res = test_halfspace_contains
% test_halfspace_contains - unit test function of containss
%
% Syntax:
%    res = test_halfspace_contains
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
% Written:       02-September-2019
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true(0);

% 2D zonotope
Z = zonotope([zeros(2,1),[1 1; -1 1]]);

% instantiate halfspaces
h_above = halfspace([0 1],3);
h_upperboundary = halfspace([0 1],2);
h_through = halfspace([0 1],1);
h_lowerboundary = halfspace([0 1],-2);
h_below = halfspace([0 1],-3);

% check if correct results for containment
% Z fully contained
res(end+1,1) = contains(h_above,Z);
% Z touching halfspace, all in
res(end+1,1) = contains(h_upperboundary,Z);
% Z partly contained
res(end+1,1) = ~contains(h_through,Z);
% Z touching halfspace, all out
res(end+1,1) = ~contains(h_lowerboundary,Z);
% Z fully outside
res(end+1,1) = ~contains(h_below,Z);


% 2D interval
I = interval([-2; -1],[1; 3]);

% instantiate halfspaces
h_above = halfspace([1 0],2);
h_upperboundary = halfspace([1 0],1);
h_through = halfspace([1 0],-1);
h_lowerboundary = halfspace([1 0],-2);
h_below = halfspace([1 0],-4);

% check if correct results for containment
% I fully contained
res(end+1,1) = contains(h_above,I);
% I touching halfspace, all in
res(end+1,1) = contains(h_upperboundary,I);
% I partly contained
res(end+1,1) = ~contains(h_through,I);
% I touching halfspace, all out
res(end+1,1) = ~contains(h_lowerboundary,I);
% I fully outside
res(end+1,1) = ~contains(h_below,I);


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
