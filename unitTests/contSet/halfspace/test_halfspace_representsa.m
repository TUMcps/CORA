function res = test_halfspace_representsa
% test_halfspace_representsa - unit test function of representsa
%
% Syntax:
%    res = test_halfspace_representsa
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
% Written:       17-September-2019
% Last update:   03-May-2020 (add empty case)
% Last revision: 20-July-2023 (MW, rename '...representsa')

% ------------------------------ BEGIN CODE -------------------------------

res = true(0);

% 1. comparison to empty set
hs = halfspace.empty(2);
res(end+1,1) = representsa(hs,'emptySet');
hs = halfspace([2;3;-1],3);
res(end+1,1) = ~representsa(hs,'emptySet');


% 2. comparison to zonotope
hs = halfspace.empty(2);
res(end+1,1) = representsa(hs,'zonotope');
hs = halfspace([2;3;-1],3);
res(end+1,1) = ~representsa(hs,'zonotope');


% 3. comparison to interval
hs = halfspace([1 0],1);
[res(end+1,1),I] = representsa(hs,'interval');
res(end+1,1) = isequal(I,interval([-Inf;-Inf],[1;Inf]));


% 4. comparison to level set
hs = halfspace([1 -1],-1);
[res(end+1,1),ls] = representsa(hs,'levelSet');
syms x y; eq = x - y + 1;
ls_ = levelSet(eq,[x;y],'<=');
res(end+1,1) = isequal(ls,ls_);


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
