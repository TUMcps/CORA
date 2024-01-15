function res = test_conHyperplane_isIntersecting
% test_conHyperplane_isIntersecting - unit test function of isIntersecting
%
% Syntax:
%    res = test_conHyperplane_isIntersecting
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
% Last revision: 09-January-2024 (MW, reformat)

% ------------------------------ BEGIN CODE -------------------------------

res = true(0);

% 2D zonotope
Z = zonotope([0;0],[1 1; -1 1]);

% instantiate halfspaces
hyp_above = conHyperplane([0 1],3);
hyp_upperboundary = conHyperplane([0 1],2);
hyp_through = conHyperplane([0 1],1);
hyp_lowerboundary = conHyperplane([0 1],-2);
hyp_below = conHyperplane([0 1],-3);

% Z fully contained
res(end+1,1) = ~isIntersecting(hyp_above,Z);
% Z touching halfspace on upper boundary (all contained)
res(end+1,1) = isIntersecting(hyp_upperboundary,Z);
% Z partly contained, intersected
res(end+1,1) = isIntersecting(hyp_through,Z);
% Z touching halfspace on lower boundary (all out)
res(end+1,1) = isIntersecting(hyp_lowerboundary,Z);
% Z fully outside
res(end+1,1) = ~isIntersecting(hyp_below,Z);


% 2D interval
I = interval([-2; -1],[1; 3]);

% instantiate halfspaces
hyp_above = conHyperplane([1 0],2);
hyp_upperboundary = conHyperplane([1 0],1);
hyp_through = conHyperplane([1 0],-1);
hyp_lowerboundary = conHyperplane([1 0],-2);
hyp_below = conHyperplane([1 0],-4);

% I fully contained
res(end+1,1) = ~isIntersecting(hyp_above,I);
% I touching halfspace on upper boundary (all contained)
res(end+1,1) = isIntersecting(hyp_upperboundary,I);
% I partly contained, intersected
res(end+1,1) = isIntersecting(hyp_through,I);
% I touching halfspace on lower boundary (all out)
res(end+1,1) = isIntersecting(hyp_lowerboundary,I);
% I fully outside
res(end+1,1) = ~isIntersecting(hyp_below,I);


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
