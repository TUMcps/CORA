function res = test_halfspace_isIntersecting
% test_halfspace_isIntersecting - unit test function of isIntersecting
%
% Syntax:
%    res = test_halfspace_isIntersecting
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
% Last update:   04-May-2020 (adapt acc. to new definition of isIntersecting)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true(0);

% intersection with empty set
% res_e = ~isIntersecting(halfspace([1 1],2),zonotope.empty(2));

% 2D zonotope
Z = zonotope([zeros(2,1),[1 1; -1 1]]);

% instantiate halfspaces
h_above = halfspace([0 1],3);
h_upperboundary = halfspace([0 1],2);
h_through = halfspace([0 1],1);
h_lowerboundary = halfspace([0 1],-2);
h_below = halfspace([0 1],-3);

% check if correct results for intersection
% Z fully contained
res(end+1,1) = isIntersecting(h_above,Z);
% Z touching halfspace on upper boundary (all contained)
res(end+1,1) = isIntersecting(h_upperboundary,Z);
% Z partly contained, intersected
res(end+1,1) = isIntersecting(h_through,Z);
% Z touching halfspace on lower boundary (all out)
res(end+1,1) = isIntersecting(h_lowerboundary,Z);
% Z fully outside
res(end+1,1) = ~isIntersecting(h_below,Z);


% 2D interval
I = interval([-2; -1],[1; 3]);

% instantiate halfspaces
h_above = halfspace([1 0],2);
h_upperboundary = halfspace([1 0],1);
h_through = halfspace([1 0],-1);
h_lowerboundary = halfspace([1 0],-2);
h_below = halfspace([1 0],-4);

% check if correct results for intersection
% I fully contained
res(end+1,1) = isIntersecting(h_above,I);
% I touching halfspace on upper boundary (all contained)
res(end+1,1) = isIntersecting(h_upperboundary,I);
% I partly contained, intersected
res(end+1,1) = isIntersecting(h_through,I);
% I touching halfspace on lower boundary (all out)
res(end+1,1) = isIntersecting(h_lowerboundary,I);
% I fully outside
res(end+1,1) = ~isIntersecting(h_below,I);


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
