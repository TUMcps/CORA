function res = test_polytope_normalizeConstraints
% test_polytope_normalizeConstraints - unit test function of normalization
%    of constraint vectors and offset
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
%                19-December-2023 (MW, ensure input polytope remains as is)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true(0);

% 1D, vertex instantiation
V = [-1, 0, 1, 2, 3];
P = polytope(V);
P_ = normalizeConstraints(P,'b');
res(end+1,1) = all(withinTol(P_.b,1) | withinTol(P_.b,0) | withinTol(P_.b,-1));
P_ = normalizeConstraints(P,'A');
res(end+1,1) = all(withinTol(vecnorm(P_.A',2,1),1));


% 2D, bounded
A = [2 1; -2 3; -2 0; -1 -4; 2 -3; 4 -1]; b = [2 3 -1 2 0 1]';
P = polytope(A,b);
% normalize offset to -1, 0, 1
P_ = normalizeConstraints(P);
res(end+1,1) = all(withinTol(P.b,b));
res(end+1,1) = all(withinTol(P_.b,1) | withinTol(P_.b,0) | withinTol(P_.b,-1));
% normalize constraints to norm 1
P_ = normalizeConstraints(P,'A');
res(end+1,1) = compareMatrices(P.A,A,0,'equal',true);
res(end+1,1) = all(withinTol(vecnorm(P_.A',2),1));


% 2D, equality constraints
Ae = [2 1; -2 3; -2 0; -1 -4; 2 -3; 4 -1]; be = [2 3 -1 2 0 1]';
P = polytope([],[],Ae,be);
% normalize offset in constraints to 0, 1
P_ = normalizeConstraints(P);
res(end+1,1) = all(withinTol(P.be,be));
res(end+1,1) = all( withinTol(P_.be,1) | withinTol(P_.be,0) );
% normalize constraints in constraints to length 1
P_ = normalizeConstraints(P,'A');
res(end+1,1) = compareMatrices(P.Ae,Ae,0,'equal',true);
res(end+1,1) = all(withinTol(vecnorm(P_.Ae',2),1));


% 2D, empty equality constraint
Ae = [0 0]; be = -1;
P = polytope([],[],Ae,be);
% normalize offset in constraints to 0 and 1
P_ = normalizeConstraints(P);
res(end+1,1) = all(withinTol(P.Ae,Ae)) && all(withinTol(P.be,be));
res(end+1,1) = all(withinTol(P_.Ae,Ae)) && all(withinTol(P_.be,1));


% 2D, degenerate using equality constraint
A = [0 1; 0 -1]; b = [3;-1]; Ae = [1 0]; be = 0;
P = polytope(A,b,Ae,be);
% normalize offset vector in constraints to -1,0,1
P_ = normalizeConstraints(P,'b');
res(end+1,1) = all(withinTol(P.b,b)) && all(withinTol(P.be,be));
res(end+1,1) = all(withinTol(P_.b,1) | withinTol(P_.b,0) | withinTol(P_.b,-1)) ...
    && all(withinTol(P_.be,1) | withinTol(P_.be,0));

% 2D, vertex instantiation
V = [1 3; 2 4; -1 -2; -1 0; 3 -2; 4 1; -3 -1]';
P = polytope(V);
P_ = normalizeConstraints(P,'b');
res(end+1,1) = all(withinTol(P_.b,1) | withinTol(P_.b,0) | withinTol(P_.b,-1));
P_ = normalizeConstraints(P,'A');
res(end+1,1) = all(withinTol(vecnorm(P_.A',2,1),1));


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
