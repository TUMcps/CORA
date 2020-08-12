function res = test_halfspace_in
% test_halfspace_in - unit test function of in
%
% Syntax:  
%    res = test_halfspace_in
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
h_above = halfspace([0 1],3);
h_upperboundary = halfspace([0 1],2);
h_through = halfspace([0 1],1);
h_lowerboundary = halfspace([0 1],-2);
h_below = halfspace([0 1],-3);

% check if correct results for containment
% Z fully contained
res_above = in(h_above,Z);
% Z touching halfspace, all in
res_upperboundary = in(h_upperboundary,Z);
% Z partly contained
res_through = in(h_through,Z);
% Z touching halfspace, all out
res_lowerboundary = in(h_lowerboundary,Z);
% Z fully outside
res_below = in(h_below,Z);

% compare results
res_zon = res_above && res_upperboundary && ~res_through && ...
    ~res_lowerboundary && ~res_below;
% -------------------------------------------------------------------------

% TEST 2: Intervals -------------------------------------------------------
% instantiate zonotope
I = interval([-2; -1],[1; 3]);
% plot(I);

% instantiate halfspaces
h_above = halfspace([1 0],2);
h_upperboundary = halfspace([1 0],1);
h_through = halfspace([1 0],-1);
h_lowerboundary = halfspace([1 0],-2);
h_below = halfspace([1 0],-4);

% check if correct results for containment
% I fully contained
res_above = in(h_above,I);
% I touching halfspace, all in
res_upperboundary = in(h_upperboundary,I);
% I partly contained
res_through = in(h_through,I);
% I touching halfspace, all out
res_lowerboundary = in(h_lowerboundary,I);
% I fully outside
res_below = in(h_below,I);

% compare results
res_int = res_above && res_upperboundary && ~res_through && ...
    ~res_lowerboundary && ~res_below;
% -------------------------------------------------------------------------


% combine tests
res = res_zon && res_int;

if res
    disp('test_in successful');
else
    disp('test_in failed');
end

%------------- END OF CODE --------------