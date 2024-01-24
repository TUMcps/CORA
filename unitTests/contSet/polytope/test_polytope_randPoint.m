function res = test_polytope_randPoint
% test_polytope_randPoint - unit test function of randPoint
%
% Syntax:
%    res = test_polytope_randPoint
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
% Written:       04-December-2022
% Last update:   03-January-2024 (MW, more tests)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true(0);
tol = 1e-6;

% empty set
P = polytope.empty(2);
p = randPoint(P);
res(end+1,1) = isnumeric(p) && isempty(p) && all(size(p) == [2,0]);

% 1D, bounded
A = [1;-1]; b = [4;-2];
P = polytope(A,b);
p = randPoint(P);
p(:,end+1) = randPoint(P,1,'extreme');
res(end+1,1) = all(contains(P,p,'exact',tol));

% 1D, single point
Ae = 2; be = 4;
P = polytope([],[],Ae,be);
p = randPoint(P);
res(end+1,1) = withinTol(p,2);

% 1D, unbounded
% A = 1; b = 2;
% P = polytope(A,b);
% p = randPoint(P);
% res(end+1,1) = p <= 2;

% 1D, fully empty
A = zeros(0,1); b = zeros(0,0);
P = polytope(A,b);
p = randPoint(P,10);
res(end+1,1) = all(contains(P,p,'exact',tol));


% 2D, bounded
A = [2 1; -1 2; -2 -3; 3 -1]; b = ones(4,1);
P = polytope(A,b);
% sample random point with different syntaxes
p = randPoint(P);
p(:,end+1) = randPoint(P,1);
p(:,end+1) = randPoint(P,1,'standard');
p(:,end+1) = randPoint(P,1,'extreme');
p(:,end+1:end+5) = randPoint(P,5,'standard');
p(:,end+1:end+5) = randPoint(P,5,'extreme');
% all have to be contained in P
res(end+1,1) = all(contains(P,p,'exact',tol));

% 2D, bounded, degenerate
A = [1 0; -1 0; 0 1; 0 -1]; b = [1;1;1;-1];
P = polytope(A,b);
p = randPoint(P);
res(end+1,1) = contains(P,p,'exact',tol);

% 2D, unbounded
A = [1 0; -1 0; 0 1]; b = [1;1;1];
P = polytope(A,b);
p = randPoint(P);
res(end+1,1) = contains(P,p);
p = randPoint(P,1,'extreme');
res(end+1,1) = contains(P,p);


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
