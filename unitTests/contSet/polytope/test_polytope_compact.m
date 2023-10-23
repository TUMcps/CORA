function res = test_polytope_compact
% test_polytope_compact - unit test function of compact
%
% Syntax:
%    res = test_polytope_compact
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

% Authors:       Mark Wetzlinger, Viktor Kotsev
% Written:       04-December-2022
% Last update:   25-May-2022 (MW, add more cases)
% Last revision: 31-July-2023 (MW, rename '...compact', merge with minVRep)

% ------------------------------ BEGIN CODE -------------------------------

res = [];

% Empty case --------------------------------------------------------------
P = polytope([1 0; -1 0],[1; -2]);

% compute minimal H-representation
P_ = compact(P);

% should return an empty object
res(end+1,1) = isemptyobject(P_);

% equality constraints
P = polytope([],[],[1 0; 0 1; 0 1],[1; -1; 0]);

% compute minimal H-representation
P_ = compact(P);

% should return an empty object
res(end+1,1) = isemptyobject(P_);


% 1D cases ----------------------------------------------------------------
% unbounded toward -Inf
P = polytope([1; 1; 1],[1;2;3]);
P_ = compact(P);
res(end+1,1) = P == P_ && length(P_.b) == 1;

% bounded
P = polytope([1; -1; 1],[1;2;3]);
P_ = compact(P);
res(end+1,1) = P == P_ && length(P_.b) == 2;

% unbounded toward Inf
P = polytope([-1; -1; -1],[1;2;3]);
P_ = compact(P);
res(end+1,1) = P == P_ && length(P_.b) == 1;

% x <= {4,5,6}, x >= {1,2,3}
% -->  x <= 4 && -x <= -3
A = [1; 1; 1; -1; -1; -1];
b = [6; 5; 4; -1; -2; -3];
P = polytope(A,b);
P_ = compact(P);

% compare to true result
A_ = [1; -1]; b_ = [4; -3];
res(end+1,1) = compareMatrices([P_.A'; P_.b'],[A_'; b_']);

% set containing constraints x <= 0.25, x >= 0.5 -> infeasible
A = [2; 0.5; 4; -2; -0.5; -1];
b = [3; 2; 1; -1; -2; -3];
P = polytope(A,b);
P_ = compact(P);

% compare results
res(end+1,1) = isemptyobject(P_);

% equality constraints x = 2, x = 1 -> infeasible
Ae = [1; 1];
be = [2; 1];
P = polytope([],[],Ae,be);
P_ = compact(P);

% compare results
res(end+1,1) = isemptyobject(P_);


% inequality constraints x <= {4,5}, x >= {0,1}
% and equality constraint x = 2 
A = [1; 1; -1; -1];
b = [4; 5; 0; -1];
Ae = 1;
be = 2;
P = polytope(A,b,Ae,be);
P_ = compact(P);

% compare results
res(end+1,1) = isempty(P_.A) && length(P.Ae) == 1;


% irredundant vertex representation
V = [1,2,3,4,5];
P = polytope(V);
P_ = compact(P,'V');

% expected result
P_true = polytope([1,5]);

% compare vertice matrices
res(end+1,1) = compareMatrices(P_.V.val, P_true.V.val);


% 2D cases ----------------------------------------------------------------

% single halfspace
P = polytope([1 0],1);

% compute minimal H-representation
P_ = compact(P);
res(end+1,1) = length(P.b) == length(P_.b);


% only two halfspaces (unbounded)
P = polytope([1 1; -1 0],[1;1]);

% compute minimal H-representation
P_ = compact(P);
res(end+1,1) = length(P.b) == length(P_.b);


% two halfspaces (parallel, feasible)
P = polytope([-3 4; -6 8],[2; 3]);

% compute minimal H-representation
P_ = compact(P);
res(end+1,1) = length(P_.b) == 1 ...
    && all(withinTol(P_.A,[-3/5 4/5])) && withinTol(P_.b,3/10);


% two halfspaces (anti-parallel, feasible)
P = polytope([-3 4; 6 -8],[2; 3]);

% compute minimal H-representation
P_ = compact(P);
res(end+1,1) = length(P_.b) == 2;

% two halfspaces (anti-parallel, infeasible)
P = polytope([-3 4; 6 -8],[-2; 3]);

% compute minimal H-representation
P_ = compact(P);
res(end+1,1) = isemptyobject(P_);


% with parallel constraints
P = polytope([1 0; -1 0; 1 0; 1 0; 1 -1; 1 -1],[1;2;3;4;5;6]);

% compute minimal H-representation
P_ = compact(P);
res(end+1,1) = length(P_.b) == 3;

% polytope without redundancies
A = [2 1; -1 3; -2 -2; 1 -4];
b = ones(4,1);
P = polytope(A,b);

% compute minimal H-representation
P_ = compact(P);
res(end+1,1) = length(P.b) == length(P_.b);


% polytope with one redundant halfspace
A = [2 1; -1 3; -2 -2; 1 -4; 1 0];
b = ones(5,1);
Pred = polytope(A,b);

% compute minimal H-representation
P_ = compact(Pred);
res(end+1,1) = size(P_.A,1) == 4 && length(P_.b) == 4;


% with equality constraints
A = [-1 0; 0 1]; b = zeros(2,1);
Ae = [1 0]; be = 0;
P = polytope(A,b,Ae,be);
P_ = compact(P);
res(end+1,1) = isequal(P,P_);


% polytope as two boxes, one fully containing the other
n = 5;
A = [eye(n);-eye(n);eye(n);-eye(n)];
b = [ones(2*n,1);2*ones(2*n,1)];
P = polytope(A,b);

% compute minimal H-representation
P_ = compact(P);
res(end+1,1) = size(P_.A,1) == 2*n && length(P_.b) == 2*n;


% more complex polytope with redundant halfspaces
A = [-1 1; 0 1; 1 0; -1 -2; -1/sqrt(2) -1/sqrt(2)];
b = [0; 1; 2; 0; -sqrt(2)];
P = polytope(A,b);

% compute minimal H-representation
P_ = compact(P);

% true solution
Pmin = polytope([-1/sqrt(2) -1/sqrt(2); 0 1; 1 0],[-sqrt(2); 1; 2]);
res(end+1,1) = P_ == Pmin;

% visualization
% figure; hold on;
% plot(P);
% plot(Pmin,[1,2],'g:');
% plot(P_,[1,2],'r--');
% close;

% unbounded degenerate polytope
A = [1 0; -1 0; 0 1; 2 0; -2 0; 0 2];
b = [1; -1; 5; 10; -0.3; 12];
P = polytope(A,b);

% compute minimal H-representation
P_ = compact(P);

% true solution
Pmin = polytope([1 0; -1 0; 0 1],[1;-1; 5]);
res(end+1,1) = P_ == Pmin;


% normal vectors only in one half of the plane
A = [1 0; 0.1 0.1; 0.1 0.2; 0.1 0.4; 0 1; -0.1 0.2; -0.2 0.1; -0.5 0.1; -1 0];
b = ones(9,1);
P = polytope(A,b);

% compute minimal H-representation
P_ = compact(P);

% only axis-aligned halfspaces should remain
res(end+1,1) = compareMatrices([1 0; -1 0; 0 1]',P_.A',1e-14);


% irredundant vertex representation
V = [2 0; 0 0; 1 0; 0 1; 2 1; 2 0; 1.5 0]';
P = polytope(V);
P_ = compact(P,'V');
% compare to true result
P_true = polytope([2 0; 0 0; 0 1; 2 1]');
res(end+1,1) = P_ == P_true;


% 3D cases ----------------------------------------------------------------

% box (no linprog in check)
A = [eye(3); -eye(3)];
b = ones(6,1);
P = polytope(A,b);

% compute minimal H-representation
P_ = compact(P);
res(end+1,1) = P == P_;

% sequence of functions
% 2D
A = [1 1; -2 1; -4 -2; 2 -3];
b = ones(4,1);
P = polytope(A,b);
P = compact(P);
res(end+1,1) = ~isempty(P.minHRep.val) && P.minHRep.val;
res(end+1,1) = isempty(P.minVRep.val);

% Boundry box should be in minimal H-representation
B = box(P);
res(end+1,1) = ~isempty(B.minHRep.val) && B.minHRep.val;

% After computing vertices, V-rep should be unknown if they were not present before
vertices(P);
res(end+1,1) = isempty(P.minVRep.val);

P = compact(P, 'V');
res(end+1,1) = ~isempty(P.minVRep.val) && P.minVRep.val;

% Now that we have the vertices the property should not change after computation
vertices(P);
res(end+1,1) = ~isempty(P.minHRep.val) && P.minHRep.val;
res(end+1,1) = ~isempty(P.minVRep.val) && P.minVRep.val;

% Intersection results in unknown properties
P2 = polytope([1 0; -1 0; 0 1; 0 -1],[1;1;1;1]);
P = P & P2;
res(end+1,1) = isempty(P.minHRep.val);
res(end+1,1) = isempty(P.minVRep.val);

% Computing minimal representation on empty V-set does not result in minimal V-representation(unless P is empty)
P = compact(P, 'all');
P = compact(P, 'V');
res(end+1,1) = ~isempty(P.minHRep.val) && P.minHRep.val;
res(end+1,1) = isempty(P.minVRep.val);

% volume implicitely computes vertices
volume(P);
P = compact(P, 'V');
res(end+1,1) = ~isempty(P.minVRep.val) && P.minVRep.val;

% No change
P = normalizeConstraints(P);
res(end+1,1) = ~isempty(P.minHRep.val) && P.minHRep.val;
res(end+1,1) = ~isempty(P.minVRep.val) && P.minVRep.val;

% Minkowski sum implicitly computes minimal H-representation (but not min V-rep)
P3 = polytope.generateRandom('Dimension',2);
P = P + P3;
res(end+1,1) = ~isempty(P.minHRep.val) && P.minHRep.val;
res(end+1,1) = isempty(P.minVRep.val);

% Minkowski difference does not (results in unknown properties)
P = minkDiff(P,P2);
res(end+1,1) = isempty(P.minHRep.val);
res(end+1,1) = isempty(P.minVRep.val);

% Project implicitly computes min H-rep
P4 = polytope.generateRandom('Dimension', 3);
P = project(P4,[2,3]);
res(end+1,1) = ~isempty(P.minHRep.val) && P.minHRep.val;
res(end+1,1) = isempty(P.minVRep.val);

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
