function res = test_zonotope_and
% test_zonotope_and - unit test function of and
%
% Syntax:  
%    res = test_zonotope_and
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

% Author:       Mark Wetzlinger
% Written:      09-September-2020
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

res = true;

% 1D (convertible to intervals)
Z1 = zonotope(0,3);
Z2 = zonotope(5,1);
Z1and2 = Z1 & Z2; % empty
res(1) = isempty(Z1and2);

Z1 = zonotope(0,4);
Z1and2 = Z1 & Z2; % only one point
res(2) = isequal(Z1and2,zonotope(4));

Z1 = zonotope(0,5);
Z1and2 = Z1 & Z2; % full-dimensional zonotope
res(3) = isequal(Z1and2,zonotope(4.5,0.5));

% 2D
Z1 = zonotope(zeros(2,1),rand(2,2));
Z2 = zonotope(5*ones(2,1),rand(2,2));
Z1and2 = Z1 & Z2; % empty
res(4) = isempty(Z1and2);

Z1 = zonotope([0;0],[1 0.5; 0 1]);
Z2 = zonotope([2.5;2.5],[1 0; 0.5 1]);
Z1and2 = Z1 & Z2; % only one point (exactly)
res(5) = ~isempty(Z1and2);

Z1 = zonotope(zeros(2,1),rand(2,2));
Z2 = zonotope(zeros(2,1),rand(2,2));
Z1and2 = Z1 & Z2; % full-dimensional intersection
res(6) = ~isempty(Z1and2);

% empty set
% Z_e = zonotope();
% res(7) = isempty(Z1 & Z_e);

res = all(res);

%------------- END OF CODE --------------
