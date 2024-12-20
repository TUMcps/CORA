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

% 1D cases ----------------------------------------------------------------

% 1D, no constraints
A = zeros(0,1); b = zeros(0,0);
P = polytope(A,b);
P_min = compact(P);
assert(representsa(P_min,'fullspace'));

% 1D, one trivially fulfilled constraint
P = polytope.Inf(1);
P_min = compact(P);
assert(representsa(P_min,'fullspace'));

% 1D, unbounded toward -Inf
A = [1; 1; 1]; b = [1;2;3];
P = polytope(A,b);
P_min = compact(P);
assert(P == P_min && isscalar(P_min.b));

% 1D, unbounded toward Inf
A = [-1; -1; -1]; b = [1;2;3];
P = polytope(A,b);
P_min = compact(P);
assert(P == P_min && isscalar(P_min.b));

% 1D, bounded, one redundant halfspace
A = [1; -1; 1]; b = [1;2;3];
P = polytope(A,b);
P_min = compact(P);
assert(P == P_min && length(P_min.b) == 2);

% 1D, bounded, redundant halfspaces
A = [1; 1; 1; -1; -1; -1];
b = [6; 5; 4; -1; -2; -3];
P = polytope(A,b);
P_min = compact(P);
% compare to true constraints
A_ = [1; -1]; b_ = [4; -3];
assert(compareMatrices([P_min.A'; P_min.b'],[A_'; b_']));

% 1D, empty (inequality constraints)
A = [2; 0.5; 4; -2; -0.5; -1];
b = [3; 2; 1; -1; -2; -3];
P = polytope(A,b);
P_min = compact(P);
assert(representsa(P_min,'emptySet'));

% 1D, empty (equality constraints)
Ae = [1; 1]; be = [2; 1];
P = polytope([],[],Ae,be);
P_min = compact(P);
assert(representsa(P_min,'emptySet'));

% 1D, degenerate, redundant inequality constraints
A = [1; 1; -1; -1]; b = [4; 5; 0; -1]; Ae = 1; be = 2;
P = polytope(A,b,Ae,be);
P_min = compact(P);
assert(isempty(P_min.A) && isscalar(P.Ae));

% 1D, vertex representation
V = [1,2,3,4,5];
P = polytope(V);
P_min = compact(P,'V');
P_true = polytope([1,5]);
assert(compareMatrices(P_min.V, P_true.V));


% 2D cases ----------------------------------------------------------------

% 2D, empty
A = [1 0; -1 0]; b = [1; -2];
P = polytope(A,b);
P_min = compact(P);
assert(representsa(P_min,'emptySet'));

% 2D, empty, only equality
Ae = [1 0; 0 1; 0 1]; be = [1; -1; 0];
P = polytope([],[],Ae,be);
P_min = compact(P);
assert(representsa(P_min,'emptySet'));

% 2D, trivially fulfilled constraints
A = [0 0; 0 0]; b = [1;2]; Ae = [0 0]; be = 0;
P = polytope(A,b,Ae,be);
P_min = compact(P);
assert(representsa(P_min,'fullspace'));

% 2D, unbounded, single halfspace
A = [1 0]; b = 1;
P = polytope(A,b);
P_min = compact(P);
assert(length(P.b) == length(P_min.b));

% 2D, unbounded
A = [1 1; -1 0]; b = [1;1];
P = polytope(A,b);
P_min = compact(P);
assert(length(P.b) == length(P_min.b));

% 2D, parallel halfspaces
A = [-3 4; -6 8]; b = [2; 3];
P = polytope(A,b);
P_min = compact(P);
assert(isscalar(P_min.b) ...
    && all(withinTol(P_min.A,[-3/5 4/5])) && withinTol(P_min.b,3/10));

% 2D, anti-parallel halfspaces
A = [-3 4; 6 -8]; b = [2; 3];
P = polytope(A,b);
P_min = compact(P);
assert(length(P_min.b) == 2);

% 2D, anti-parallel halfspaces, empty
A = [-3 4; 6 -8]; b = [-2; 3];
P = polytope(A,b);
P_min = compact(P);
assert(representsa(P_min,'emptySet'));

% 2D, some parallel constraints
A = [1 0; -1 0; 1 0; 1 0; 1 -1; 1 -1]; b = [1;2;3;4;5;6];
P = polytope(A,b);
P_min = compact(P);
assert(length(P_min.b) == 3);

% 2D, no redundancies
A = [2 1; -1 3; -2 -2; 1 -4]; b = ones(4,1);
P = polytope(A,b);
P_min = compact(P);
assert(length(P.b) == length(P_min.b));

% 2D, one redundant halfspace
A = [2 1; -1 3; -2 -2; 1 -4; 1 0]; b = ones(5,1);
P = polytope(A,b);
P_min = compact(P);
assert(size(P_min.A,1) == 4 && length(P_min.b) == 4);

% 2D, with equality constraints
A = [-1 0; 0 1]; b = zeros(2,1);
Ae = [1 0]; be = 0;
P = polytope(A,b,Ae,be);
P_min = compact(P);
assert(isequal(P,P_min));

% 2D, with redundant halfspaces
A = [-1 1; 0 1; 1 0; -1 -2; -1/sqrt(2) -1/sqrt(2)];
b = [0; 1; 2; 0; -sqrt(2)];
P = polytope(A,b);
P_min = compact(P);
P_true = polytope([-1/sqrt(2) -1/sqrt(2); 0 1; 1 0],[-sqrt(2); 1; 2]);
assert(P_min == P_true);

% 2D, unbounded, degenerate
A = [1 0; -1 0; 0 1; 2 0; -2 0; 0 2];
b = [1; -1; 5; 10; -0.3; 12];
P = polytope(A,b);
P_min = compact(P);
P_true = polytope([1 0; -1 0; 0 1],[1;-1; 5]);
assert(P_min == P_true);

% 2D, normal vectors only in one half of the plane
A = [1 0; 0.1 0.1; 0.1 0.2; 0.1 0.4; 0 1; -0.1 0.2; -0.2 0.1; -0.5 0.1; -1 0];
b = ones(9,1);
P = polytope(A,b);
P_min = compact(P);
% only axis-aligned halfspaces should remain
assert(compareMatrices([1 0; -1 0; 0 1]',P_min.A',1e-14));

% 2D, vertex representation
V = [2 0; 0 0; 1 0; 0 1; 2 1; 2 0; 1.5 0]';
P = polytope(V);
P_min = compact(P,'V');
% compare to true result
V_true = [2 0; 0 0; 0 1; 2 1]';
P_true = polytope(V_true);
assert(P_min == P_true);

% 2D, one redundant constraint
A = [1 0; 0 1; 0 1]; b = [1; -1; 0];
P = polytope(A,b);
P_min = compact(P);
A_true = [1 0; 0 1]; b_true = [1; -1];
P_true = polytope(A_true,b_true);
assert(P_min == P_true);

% 2D, rewrite inequality pairs as equality constraints
A = [1 0; -1 0]; b = [1; -1];
P = polytope(A,b);
P_min = compact(P,'AtoAe');
assert(isempty(P_min.A) && isempty(P_min.b) ...
    && size(P_min.Ae,1) == 1 && isscalar(P_min.be) && isequal(P_min,P));


% nD cases ----------------------------------------------------------------

% 3D, box (no linprog in check), no redundancies
A = [eye(3); -eye(3)];
b = ones(6,1);
P = polytope(A,b);
P_min = compact(P);
assert(P == P_min);

% 5D, polytope as two boxes, one fully containing the other
n = 5;
A = [eye(n);-eye(n);eye(n);-eye(n)]; b = [ones(2*n,1);2*ones(2*n,1)];
P = polytope(A,b);
P_min = compact(P);
assert(size(P_min.A,1) == 2*n && length(P_min.b) == 2*n);


% sequence of functions ---------------------------------------------------
% 2D, bounded
A = [1 1; -2 1; -4 -2; 2 -3]; b = ones(4,1);
P = polytope(A,b);
P = compact(P);
% check minimal H representation
assert(~isempty(P.minHRep.val) && P.minHRep.val);
assert(isempty(P.minVRep.val));

% boundary box should be in minimal H-representation
B = box(P);
assert(~isempty(B.minHRep.val) && B.minHRep.val);

% compute vertices 
vertices(P);
P = compact(P, 'V');
assert(~isempty(P.minVRep.val) && P.minVRep.val);

% Now that we have the vertices the property should not change after
% computation
vertices(P);
assert(~isempty(P.minHRep.val) && P.minHRep.val);
assert(~isempty(P.minVRep.val) && P.minVRep.val);

% Intersection results in unknown properties
P2 = polytope([1 0; -1 0; 0 1; 0 -1],[1;1;1;1]);
P = P & P2;
assert(isempty(P.minHRep.val));
assert(isempty(P.minVRep.val));

% Computing minimal representation on empty V-set does not result in 
% minimal V-representation (unless P is empty)
P = compact(P, 'all');
vertices(P);
P = compact(P, 'V');
assert(~isempty(P.minHRep.val) && P.minHRep.val);

% volume implicitly computes vertices
volume(P);
P = compact(P, 'V');
assert(~isempty(P.minVRep.val) && P.minVRep.val);

% normalization of constraints does not cause any changes
P = normalizeConstraints(P);
assert(~isempty(P.minHRep.val) && P.minHRep.val);
assert(~isempty(P.minVRep.val) && P.minVRep.val);

% Minkowski sum (using project) does not implicitly compute minimal
% H-representation anymore
A = [ -0.793 -0.609 ; -0.609 0.793 ; 0.992 -0.130 ; -0.942 -0.335 ; 0.729 -0.684 ; 0.547 0.837 ];
b = [ 1.027 ; -0.314 ; 1.910 ; 0.620 ; 2.292 ; 0.541 ];
P3 = polytope(A,b);
P = P + P3;
assert(isempty(P.minHRep.val));
assert(isempty(P.minVRep.val));

% Minkowski difference does neither (results in unknown properties)
P = minkDiff(P,P2);
assert(isempty(P.minHRep.val));
assert(isempty(P.minVRep.val));

% Project does not implicitly compute min H-representation anymore
A = [ 0.194 0.615 -0.764 ; -0.757 0.589 0.282 ; 0.624 0.524 0.580 ; -0.035 -0.998 -0.056 ; -0.569 -0.716 -0.404 ; 0.923 -0.179 -0.342 ; -0.807 0.097 0.583 ];
b = [ 2.206 ; 1.319 ; 1.202 ; 0.003 ; 0.480 ; 1.169 ; 0.588 ];
P4 = polytope(A,b);
P = project(P4,[2,3]);
assert(isempty(P.minHRep.val));
assert(isempty(P.minVRep.val));

% test empty polytope
P = polytope([],[],0,1);
P = compact(P);
assert(representsa_(P,'emptySet'),1e-12);

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
