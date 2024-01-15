function res = test_polytope_polytope
% test_polytope_polytope - unit test function of polytope (constructor)
%
% Syntax:
%    res = test_polytope_polytope
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

% Authors:       Viktor Kotsev, Mark Wetzlinger
% Written:       25-April-2022
% Last update:   25-July-2023 (MW, integrate computeHRep, more tests)
%                08-December-2023 (MW, unbounded cases)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true(0);
tol = 1e-12;

% empty object
P = polytope.empty(2);
res(end+1,1) = representsa(P,'emptySet') ...
                && ~isempty(P.emptySet.val) && P.emptySet.val ...
                && ~isempty(P.bounded.val) && P.bounded.val ...
                && ~isempty(P.fullDim.val) && ~P.fullDim.val;


% instantiate H-representation --------------------------------------------

% only inequalities
A = [1 0 -1 0 1; 0 1 0 -1 1]';
b = [3; 2; 3; 2; 1];
P = polytope(A,b);
res(end+1,1) = compareMatrices([A,b]',[P.A,P.b]',tol);

% only equalities
Ae = [1 1]; be = 1;
P = polytope([],[],Ae,be);
res(end+1,1) = compareMatrices([Ae,be]',[P.Ae,P.be]',tol);

% inequalities and equalities
A = [1 0 1; -1 0 0; 0 -1 0]; b = [1; 2; 1];
Ae = [0 1 1]; be = 0;
P = polytope(A,b,Ae,be);
res(end+1,1) = compareMatrices([A,b]',[P.A,P.b]',tol);
res(end+1,1) = compareMatrices([Ae,be]',[P.Ae,P.be]',tol);

% inequalities with Inf
A = [1 0; 0 1];
b = [Inf; 3];
P = polytope(A,b);
res(end+1,1) = compareMatrices([0,1,3]',[P.A,P.b]',tol);

% inequalities with -Inf
A = [1 0 0; -1 1 0; -1 -1 0; 0 0 1; 0 0 -1];
b = [3;2;3;-Inf;1];
P = polytope(A,b);
res(end+1,1) = representsa(P,'emptySet') ...
                && ~isempty(P.emptySet.val) && P.emptySet.val ...
                && ~isempty(P.bounded.val) && P.bounded.val ...
                && ~isempty(P.fullDim.val) && ~P.fullDim.val;

% equalities with Inf
A = [1 0; -1 1; -1 -1]; b = [1;1;1];
Ae = [1 1]; be = Inf;
P = polytope(A,b,Ae,be);
res(end+1,1) = representsa(P,'emptySet') ...
                && ~isempty(P.emptySet.val) && P.emptySet.val ...
                && ~isempty(P.bounded.val) && P.bounded.val ...
                && ~isempty(P.fullDim.val) && ~P.fullDim.val;


% instantiate from V-representation ---------------------------------------

% 1D
V = [-2 1 5 4 2 2];
% compute halfspace representation
P = polytope(V);
% true solution
P_true = polytope([1;-1],[5;2]);
res(end+1,1) = P == P_true ...
    && ~isempty(P.emptySet.val) && ~P.emptySet.val ...
    && ~isempty(P.bounded.val) && P.bounded.val;

% 1D: only one point
V = -5;
P = polytope(V);
% true solution
P_true = polytope([],[],1,-5);
res(end+1,1) = P == P_true;

% 1D: only one point
V = 3;
% compute halfspace representation
P = polytope(V);
% true solution
P_true = polytope([1;-1],[3;-3]);
res(end+1,1) = P == P_true;

% 1D, with Inf
V = [-Inf, 2];
% compute halfspace representation
P = polytope(V);
% true solution
P_true = polytope(1,2);
res(end+1,1) = P == P_true ...
    && ~isempty(P.emptySet.val) && ~P.emptySet.val ...
    && ~isempty(P.bounded.val) && ~P.bounded.val ...
    && ~isempty(P.fullDim.val) && P.fullDim.val;


% 2D: non-degenerate
V = [1 -1 -1 1; 1 -1 1 -1];
% compute halfspace representation
P = polytope(V);

% true solution
P_true = polytope([1 0; 0 1;-1 0;0 -1],ones(4,1));
res(end+1,1) = P == P_true ...
    && ~isempty(P.emptySet.val) && ~P.emptySet.val ...
    && ~isempty(P.bounded.val) && P.bounded.val;

% 2D: only one point
V = [1;-2];
% compute halfspace representation
P = polytope(V);

% true solution
P_true = polytope([1 0; 0 1;-1 0;0 -1],[1;-2;-1;2]);
res(end+1,1) = P == P_true;

% 2D: degenerate #1
V = [2 2; 2 6]';
% compute halfspace representation
P = polytope(V);

% true solution (scale offset accordingly)
P_true = polytope([1 0; 0 1; -1 0; 0 -1],[2;6;-2;-2]);
res(end+1,1) = P == P_true;

% 2D: degenerate #2
V = [2 3; -5 3]';
% compute halfspace representation
P = polytope(V);

% true solution (scale offset accordingly)
P_true = polytope([1 0; 0 1; -1 0; 0 -1],[2;3;5;-3]);
res(end+1,1) = P == P_true ...
    && ~isempty(P.emptySet.val) && ~P.emptySet.val ...
    && ~isempty(P.bounded.val) && P.bounded.val;

% 2D: degenerate #2
V = [-3 4; 0 0; 3 -4]';
% compute halfspace representation
P = polytope(V);

% true solution (scale offset accordingly)
P_true = polytope([3 -4; -3 4; 4 3; -4 -3],[25;25;0;0]);
res(end+1,1) = P == P_true ...
    && ~isempty(P.emptySet.val) && ~P.emptySet.val ...
    && ~isempty(P.bounded.val) && P.bounded.val;


% 3D: degenerate set (rotated square + translation)
shift = [10;-5;3];
V = [1 0 0; 0 1 0; -1 0 0; 0 -1 0]' + shift;
% compute halfspace representation
P = polytope(V);

% another true solution (non-unique)
P_true = polytope([1 1 0; -1 1 0; -1 -1 0; 1 -1 0],ones(4,1),[0 0 1],0) + shift;
res(end+1,1) = isequal(P,P_true,1e-8) ...
    && ~isempty(P.emptySet.val) && ~P.emptySet.val ...
    && ~isempty(P.bounded.val) && P.bounded.val;

% 3D: single vertex
V = [1; 1; 1];
% compute halfspace representation
P = polytope(V);

% another true solution (non-unique)
P_true = polytope([eye(3); -eye(3)],[ones(3,1);-ones(3,1)]);
res(end+1,1) = P == P_true ...
    && ~isempty(P.emptySet.val) && ~P.emptySet.val ...
    && ~isempty(P.bounded.val) && P.bounded.val;

% 5D: 3D cube
V = vertices(interval(-ones(3,1),ones(3,1)));
V = [V; zeros(2,2^3)];
[M,~,~] = svd(randn(5));
V = M*V;
P = polytope(V);

% another true solution (non-unique)
A_true = [eye(3) zeros(3,2); -eye(3) zeros(3,2)];
b_true = ones(6,1);
Ae_true = [0 0 0 1 0; 0 0 0 0 1]; be_true = [0;0];
P_true = M * polytope(A_true,b_true,Ae_true,be_true);
res(end+1,1) = isequal(P,P_true,1e-8); % && ~P.emptySet.val && P.bounded.val; % add empty property checks!


% copy constructor (empty)
A = [1 0; -1 0]; b = [2; -4];
P = polytope(A,b);
representsa_(P,'emptySet',1e-10);
P_copy = polytope(P);
res(end+1,1) = ~isempty(P_copy.emptySet.val) && ~isempty(P_copy.bounded.val) ...
    && ~isempty(P_copy.fullDim.val);

% copy constructor (bounded)
A = [1 0; -1 1; -1 -1]; b = [2; 2; 1];
P = polytope(A,b);
isBounded(P);
P_copy = polytope(P);
res(end+1,1) = ~isempty(P_copy.bounded.val);


% combine results
res = all(res);

% wrong initializations
A = [1 0 -1 0 1; 0 1 0 -1 1]';
b = [3; 2; 3; 2; 1];
b_ = [3; 2; 3; 2;];
Vinf = [1 Inf 0; 0 -1 1];
Vnan = [1 NaN 0; 0 -1 1];

% dimension mismatch
try
    P = polytope(A,b_);
    throw(CORAerror('CORA:testFailed'));
end

% empty argument
try
    P = polytope([],b);
    throw(CORAerror('CORA:testFailed'));
end

% too many arguments
try
    P = polytope(A,b,A,b,A);
    throw(CORAerror('CORA:testFailed'));
end

% non-finite vertices
try
    P = polytope(Vinf);
    throw(CORAerror('CORA:testFailed'));
end
try
    P = polytope(Vnan);
    throw(CORAerror('CORA:testFailed'));
end

% ------------------------------ END OF CODE ------------------------------
