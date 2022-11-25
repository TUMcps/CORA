function res = test_conHyperplane_isIntersecting
% test_isIntersecting - unit test function of isIntersecting
%
% Syntax:  
%    res = test_isIntersecting
%
% Inputs:
%    -
%
% Outputs:
%    res - boolean 
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Author:       Mark Wetzlinger
% Written:      02-Sep-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% TEST 1: Zonotopes -------------------------------------------------------
% instantiate zonotope
Z = zonotope([zeros(2,1),[1 1; -1 1]]);
% plot(Z);

% instantiate halfspaces
hyp_above = conHyperplane([0 1],3);
hyp_upperboundary = conHyperplane([0 1],2);
hyp_through = conHyperplane([0 1],1);
hyp_lowerboundary = conHyperplane([0 1],-2);
hyp_below = conHyperplane([0 1],-3);

% check if correct results for intersection
% Z fully contained
res_above = isIntersecting(hyp_above,Z);
% Z touching halfspace on upper boundary (all contained)
res_upperboundary = isIntersecting(hyp_upperboundary,Z);
% Z partly contained, intersected
res_through = isIntersecting(hyp_through,Z);
% Z touching halfspace on lower boundary (all out)
res_lowerboundary = isIntersecting(hyp_lowerboundary,Z);
% Z fully outside
res_below = isIntersecting(hyp_below,Z);

% compare results
res_zon = ~res_above && res_upperboundary && res_through && ...
    res_lowerboundary && ~res_below;
% -------------------------------------------------------------------------

% TEST 2: Intervals -------------------------------------------------------
% instantiate zonotope
I = interval([-2; -1],[1; 3]);
% plot(I);

% instantiate halfspaces
hyp_above = conHyperplane([1 0],2);
hyp_upperboundary = conHyperplane([1 0],1);
hyp_through = conHyperplane([1 0],-1);
hyp_lowerboundary = conHyperplane([1 0],-2);
hyp_below = conHyperplane([1 0],-4);

% check if correct results for intersection
% I fully contained
res_above = isIntersecting(hyp_above,I);
% I touching halfspace on upper boundary (all contained)
res_upperboundary = isIntersecting(hyp_upperboundary,I);
% I partly contained, intersected
res_through = isIntersecting(hyp_through,I);
% I touching halfspace on lower boundary (all out)
res_lowerboundary = isIntersecting(hyp_lowerboundary,I);
% I fully outside
res_below = isIntersecting(hyp_below,I);

% compare results
res_int = ~res_above && res_upperboundary && res_through && ...
    res_lowerboundary && ~res_below;
% -------------------------------------------------------------------------


% combine tests
res = res_zon && res_int;

if res
    disp('test_isIntersecting successful');
else
    disp('test_isIntersecting failed');
end

%------------- END OF CODE --------------