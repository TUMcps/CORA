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
% See also: none

% Authors:       Mark Wetzlinger
% Written:       09-September-2020
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% 1D (convertible to intervals)
Z1 = zonotope(0,3);
Z2 = zonotope(5,1);
Z_and = Z1 & Z2; % empty
assert(representsa(Z_and,'emptySet'));

Z1 = zonotope(0,4);
Z_and = Z1 & Z2; % only one point
assert(isequal(Z_and,zonotope(4)));

Z1 = zonotope(0,5);
Z_and = Z1 & Z2; % full-dimensional zonotope
assert(isequal(Z_and,zonotope(4.5,0.5)));

% 2D
Z1 = zonotope(zeros(2,1),rand(2,2));
Z2 = zonotope(5*ones(2,1),rand(2,2));
Z_and = Z1 & Z2; % empty
assert(representsa(Z_and,'emptySet'));

Z1 = zonotope([0;0],[1 0.5; 0 1]);
Z2 = zonotope([2.5;2.5],[1 0; 0.5 1]);
Z_and = Z1 & Z2; % only one point (exactly)
assert(~representsa(Z_and,'emptySet'));

Z1 = zonotope(zeros(2,1),rand(2,2));
Z2 = zonotope(zeros(2,1),rand(2,2));
Z_and = Z1 & Z2; % full-dimensional intersection
assert(~representsa(Z_and,'emptySet'));

% empty set
Z_e = zonotope.empty(2);
assert(representsa(Z1 & Z_e,'emptySet'));


% averaging method
Z1 = zonotope([1;2],[1 -1 2 0; 1 4 0 1]);
Z2 = zonotope([-3;-4],[1 0 2; -1 1 2]);
Z_and = and(Z1,Z2,'averaging');
assert(~representsa(Z_and,'emptySet'));
Z_and = and_(Z1,Z2,'averaging','normGen',false);
assert(~representsa(Z_and,'emptySet'));
Z_and = and_(Z1,Z2,'averaging','normGen',true,0.8);
assert(~representsa(Z_and,'emptySet'));
Z_and = and_(Z1,Z2,'averaging','radius');
assert(~representsa(Z_and,'emptySet'));
Z_and = and_(Z1,Z2,'averaging','volume');
assert(~representsa(Z_and,'emptySet'));


% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
