function res = test_polytope_vertices
% test_polytope_vertices - unit test function of vertices
%
% Syntax:
%    res = test_polytope_vertices
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
% Written:       29-November-2022
% Last update:   24-November-2023 (MW, check errors for unbounded cases)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true(0);

% 1D, empty, only inequalities
A = [2; 1; -4]; b = [5; 0.5; -3];
P = polytope(A,b);
V = vertices(P);
res(end+1,1) = isempty(V);

% 1D, empty, inequalities and equalities
A = [1;-4]; b = [4;-2]; Ae = 5; be = 100;
P = polytope(A,b,Ae,be);
V = vertices(P);
res(end+1,1) = isempty(V);

% 1D, bounded
A = [1; -2]; b = [1.25; 1];
P = polytope(A,b);
V = vertices(P);
V_ = [-0.5, 1.25];
res(end+1,1) = all(withinTol(V,V_));

% 1D, bounded, non-minimal representation
A = [2; 1; -1; -2]; b = [1; 2; 5; 2];
P = polytope(A,b);
V = vertices(P);
V_ = [-1 0.5];
res(end+1,1) = all(withinTol(V,V_));

% 1D, unbounded
A = 1; b = 0;
P = polytope(A,b);
V = vertices(P);
V_ = [-Inf, 0];
res(end+1,1) = all(V == V_);

% 1D, unbounded
A = 3; b = 6;
P = polytope(A,b);
V = vertices(P);
V_ = [-Inf, 2];
res(end+1,1) = all(V == V_);

% 1D, degenerate, only inequalities
A = [2; -1]; b = [4; -2];
P = polytope(A,b);
V = vertices(P);
V_ = 2;
res(end+1,1) = all(withinTol(V,V_));

% 1D, degenerate, equality
Ae = 1; be = 2;
P = polytope([],[],Ae,be);
V = vertices(P);
V_ = 2;
res(end+1,1) = all(withinTol(V,V_));


% 2D, bounded
A = [2 1; 2 -2; -2 -2; -1 2]; b = [2; 2; 2; 2];
P = polytope(A,b);
V = vertices(P);
V_comb = vertices(P,'comb');
V_ = [0.4 1.2; 1 0; 0 -1; -4/3 1/3]';
res(end+1,1) = compareMatrices(V,V_) && compareMatrices(V,V_comb,1e-14);

% 2D, vertex instantiation
V_ = [3 0; 2 2; -1 3; -2 0; 0 -1]';
P = polytope(V_);
V = vertices(P);
V_comb = vertices(P,'comb');
res(end+1,1) = compareMatrices(V,V_,1e-14) && compareMatrices(V,V_comb,1e-14);

% 2D, bounded, degenerate (single point)
A = [1 1; 1 -1; -1 0]; b = zeros(3,1);
P = polytope(A,b);
V = vertices(P);
V_comb = vertices(P,'comb');
V_ = [0;0];
res(end+1,1) = compareMatrices(V,V_,1e-14) && compareMatrices(V,V_comb,1e-14);

% 2D, bounded, degenerate (line)
A = [1 1; 1 -1; -1 -1; -1 1]; b = [1; 0; 1; 0];
P = polytope(A,b);
V = vertices(P);
V_comb = vertices(P,'comb');
V_ = [0.5 0.5; -0.5 -0.5]';
res(end+1,1) = compareMatrices(V,V_,1e-14) && compareMatrices(V,V_comb,1e-14);


% 3D, degenerate (2D simplex)
A = [-1 0 0; 0 -1 0; 1 1 0]; b = [0;0;2];
Ae = [0 0 1]; be = 0;
P = polytope(A,b,Ae,be);
V = vertices(P);
V_ = [2 0 0; 0 2 0; 0 0 0]';
res(end+1,1) = compareMatrices(V,V_,1e-14);

% 3D, unit box
n = 3; A = [eye(n); -eye(n)]; b = [ones(2*n,1)];
P = polytope(A,b);
V = vertices(P);
V_ = [1 1 1; 1 1 -1; 1 -1 1; -1 1 1; 1 -1 -1; -1 1 -1; -1 -1 1; -1 -1 -1]';
res(end+1,1) = compareMatrices(V,V_);

% 3D, degenerate unit box: square 
n = 3; A = [eye(n); -eye(n)]; b = [1; 1; 0; 1; 1; 0];
P = polytope(A,b);
V = vertices(P);
V_comb = vertices(P,'comb');
V_ = [1 1 0; -1 1 0; -1 -1 0; 1 -1 0]';
res(end+1,1) = compareMatrices(V,V_,1e-14) && compareMatrices(V,V_comb,1e-14);


% 4D, degenerate (unit square)
A = [eye(2) zeros(2); -eye(2) zeros(2)]; b = [ones(4,1)];
Ae = [0 0 1 0; 0 0 0 1]; be = [0;0];
P = polytope(A,b,Ae,be);
V = vertices(P);
V_ = [1 1 0 0; 1 -1 0 0; -1 1 0 0; -1 -1 0 0]';
res(end+1,1) = compareMatrices(V,V_,1e-14);

% 4D, degenerate (rotated unit square)
A = [eye(2) zeros(2); -eye(2) zeros(2)]; b = [ones(4,1)];
Ae = [0 0 1 0; 0 0 0 1]; be = [0;0];
P = polytope(A,b,Ae,be);
M = [1 3 -2 4; 3 -2 4 -1; 3 -2 -1 3; 4 -3 -2 1];
V = vertices(M*P);
res(end+1,1) = compareMatrices(V,M*V_,1e-14);


% 7D, degenerate (subspace computation required)
A = [0       1  0   0  0  0   0;
     0      -1  0   0  0  0   0;
     0       0  0   1  0  0   0;
     0       0  0  -1  0  0   0;
     0       0  0   0  1  0   0;
     0       0  0   0 -1  0   0;
     0       0  0   0  0  0   1;
     0       0  0   0  0  0  -1;
     sqrt(2) 0 -0.5 0  0  0.5 0;
     sqrt(2) 0  0.5 0  0 -0.5 0;
     0       0  0   0  0  1   0;
     0       0 -1   0  0  0   0;
    -1       0  0   0  0  0   0];
b = [0;0;1;-1;0;0;0;0;sqrt(2);sqrt(2);0.26;-0.25;-1];
P = polytope(A,b);
V = vertices(P);
V_ = [1 0 0.25 1 0 0.25 0; 1 0 0.26 1 0 0.26 0]';
res(end+1,1) = compareMatrices(V,V_,1e-10);


% combine results
res = all(res);


% 2D, unbounded (should throw an error!)
A = [1 0]; b = 1;
P = polytope(A,b);
try
    V = vertices(P);
    res = false;
end

% 2D, unbounded (should throw an error!)
A = [1 0; -1 0]; b = [1;1];
P = polytope(A,b);
try
    V = vertices(P);
    res = false;
end
[M,~,~] = svd(rand(2));
P_ = M*P;
try
    V = vertices(P_);
    res = false;
end

% 2D, unbounded (should throw an error!)
A = [1 0; -1 0; 0 -1]; b = [1;1;1];
P = polytope(A,b);
try
    V = vertices(P);
    res = false;
end

% 2D, unbounded but not axis-aligned (should throw an error!)
A = [1 -0.1; 0.1 -1; -0.1 -1; -1 -0.1]; b = ones(4,1);
P = polytope(A,b);
try
    V = vertices(P);
    res = false;
end

% ------------------------------ END OF CODE ------------------------------
