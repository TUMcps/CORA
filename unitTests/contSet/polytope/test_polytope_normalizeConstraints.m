function res = test_polytope_normalizeConstraints
% test_polytope_normalizeConstraints - unit test function of normlization
%
% Syntax:
%    res = test_polytope_normalizeConstraints
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
% Written:       01-December-2022
% Last update:   13-December-2022 (MW, add cases for type = 'A')
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = [];

% init polytope
A = [2 1; -2 3; -2 0; -1 -4; 2 -3; 4 -1];
b = [2 3 -1 2 0 1]';
P = polytope(A,b);

% normalize offset in constraints
P_ = normalizeConstraints(P);

% check if all inequality constraints are 1
res(end+1,1) = all(withinTol(P_.b,1) | withinTol(P_.b,0) | withinTol(P_.b,-1));

% normalize constraints in constraints
P_ = normalizeConstraints(P,'A');

% check if all inequality constraints have norm 1
res(end+1,1) = all(withinTol(vecnorm(P_.A',2),1));


% init polytope with equality constraints
Ae = [2 1; -2 3; -2 0; -1 -4; 2 -3; 4 -1];
be = [2 3 -1 2 0 1]';
P = polytope([],[],Ae,be);

% normalize offset in constraints
P_ = normalizeConstraints(P);

% check normalized values
res(end+1,1) = all( withinTol(P_.be,1) | withinTol(P_.be,0) );

% normalize constraints in constraints
P_ = normalizeConstraints(P,'A');

% check if all equality constraints have norm 1
res(end+1,1) = all(withinTol(vecnorm(P_.Ae',2),1));


% 2D, degenerate using equality constraint
A = [0 1; 0 -1]; b = [3;-1];
Ae = [1 0]; be = 0;
P = polytope(A,b,Ae,be);

% normalize offset vector in constraints
P_ = normalizeConstraints(P,'b');

% check if offset vectors are -1|0|1
res(end+1,1) = all(withinTol(P_.b,1) | withinTol(P_.b,0) | withinTol(P_.b,-1)) ...
    && all(withinTol(P_.be,1) | withinTol(P_.be,0));


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
