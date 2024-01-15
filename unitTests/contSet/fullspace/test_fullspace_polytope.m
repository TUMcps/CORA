function res = test_fullspace_polytope
% test_fullspace_polytope - unit test function of polytope conversion
%
% Syntax:
%    res = test_fullspace_polytope
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
% Written:       15-December-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true(0);

% init fullspace
fs = fullspace(2);
P = polytope(fs);
res(end+1,1) = ~isBounded(P);
res(end+1,1) = supportFunc(P,[1;0],'upper') == Inf;
res(end+1,1) = supportFunc(P,[-1;0],'upper') == Inf;
res(end+1,1) = supportFunc(P,[0;1],'upper') == Inf;
res(end+1,1) = supportFunc(P,[0;-1],'upper') == Inf;


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
