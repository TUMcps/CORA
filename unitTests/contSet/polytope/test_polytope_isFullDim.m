function res = test_polytope_isFullDim
% test_polytope_isFullDim - unit test function of isFullDim
%
% Syntax:
%    res = test_polytope_isFullDim
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

% Authors:       Viktor Kotsev, Adrian Kulmburg, Mark Wetzlinger
% Written:       09-May-2022
% Last update:   25-May-2023 (AK, added tests for subspace)
%                27-July-2023 (MW, add 1D cases)
%                12-September-2023 (TL, fixed random polytopes, added empty property checks)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true(0);

% empty object
P = polytope.empty(2);
res(end+1,1) = ~isFullDim(P) && ~isempty(P.fullDim.val) && ~P.fullDim.val;

% fullspace object
P = polytope.Inf(2);
res(end+1,1) = isFullDim(P) && ~isempty(P.fullDim.val) && P.fullDim.val;


% 1D, only inequalities, bounded
A = [2;-1]; b = [6;1];
P = polytope(A,b);
res(end+1,1) = isFullDim(P) && ~isempty(P.fullDim.val) && P.fullDim.val;

% 1D, only equalities, single point
Ae = 3; be = 5;
P = polytope([],[],Ae,be);
res(end+1,1) = ~isFullDim(P) && ~isempty(P.fullDim.val) && ~P.fullDim.val;

% 1D, only inequalities, unbounded
A = [3;2;4]; b = [5;2;-3];
P = polytope(A,b);
res(end+1,1) = isFullDim(P) && ~isempty(P.fullDim.val) && P.fullDim.val;

% 1D, only inequalities, empty
Ae = [1;4]; be = [2;-5];
P = polytope([],[],Ae,be);
res(end+1,1) = ~isFullDim(P) && ~isempty(P.fullDim.val) && ~P.fullDim.val ...
    && ~isempty(P.emptySet.val) && P.emptySet.val;

% 1D, inequalities and equalities, empty
A = [1;-4]; b = [4;-2]; Ae = 5; be = 100;
P = polytope(A,b,Ae,be);
res(end+1,1) = ~isFullDim(P) && ~isempty(P.fullDim.val) && ~P.fullDim.val ...
    && ~isempty(P.emptySet.val) && P.emptySet.val;

% 1D, fully empty
A = zeros(0,1); b = zeros(0,0);
P = polytope(A,b);
res(end+1,1) = isFullDim(P) && ~isempty(P.fullDim.val) && P.fullDim.val ...
    && ~isempty(P.emptySet.val) && ~P.emptySet.val;


% 2D, empty
A = [1 0]; b = 3; Ae = [1 0]; be = 4;
P = polytope(A,b,Ae,be);
res(end+1,1) = ~isFullDim(P) && ~isempty(P.fullDim.val) && ~P.fullDim.val;

% 2D, non-degenerate, vertex instantiation
V = [2 0; -2 0; 0 2; 0 -2]';
P = polytope(V);
res(end+1,1) = isFullDim(P) && ~isempty(P.fullDim.val) && P.fullDim.val;

% 2D, non-degenerate, bounded
A = [-1 -1; 1 0;-1 0; 0 1; 0 -1]; b = [2; 3; 2; 3; 2];
P = polytope(A,b);
res(end+1,1) = isFullDim(P) && ~isempty(P.fullDim.val) && P.fullDim.val;

% 2D, degenerate
A = [1 0; -1 0; 0 1; 0 -1]; b = [2; 2; 2; -2];
P = polytope(A,b);
res(end+1,1) = ~isFullDim(P) && ~isempty(P.fullDim.val) && ~P.fullDim.val;

% 2D, degenerate
A = [1 1; 1 -1; -1 0]; b = zeros(3,1);
P = polytope(A,b);
res(end+1,1) = ~isFullDim(P) && ~isempty(P.fullDim.val) && ~P.fullDim.val;

% 2D, non-degenerate, unbounded
A = [-1 0; 0 -1]; b = [-1; -1];
P = polytope(A,b);
res(end+1,1) = isFullDim(P) && ~isempty(P.fullDim.val) && P.fullDim.val;

% 2D, degenerate, unbounded
Ae = [1 0]; be = 0;
P = polytope([],[],Ae,be);
res(end+1,1) = ~isFullDim(P) && ~isempty(P.fullDim.val) && ~P.fullDim.val;

% 2D, degenerate, bounded, inequality and equality constraints
A = [-1 0; 0 -1]; b = [-1; -1]; Ae = [1 0]; be = 0;
P = polytope(A,b,Ae,be);
res(end+1,1) = ~isFullDim(P) && ~isempty(P.fullDim.val) && ~P.fullDim.val;

% 2D, non-degenerate, unbounded (with subspace computation)
% ...unit square in x1-x2
A = [1 0;-1 0;0 1;0 -1]; b = [1;1;1;1];
P = polytope(A,b);
[res_, X] = isFullDim(P);
res(end+1,1) = res_ && rank(X) == 2 && all(size(X) == [2,2]);


% 3D, degenerate, unbounded
A = [1 1 0; 1 -1 0; -1 0 0]; b = zeros(3,1);
P = polytope(A,b);
res(end+1,1) = ~isFullDim(P) && ~isempty(P.fullDim.val) && ~P.fullDim.val;

% 3D, degenerate, bounded
A = [1 1 0; 1 -1 0; -1 0 0; 0 0 1; 0 0 -1]; b = [0; 0; 0; 1; 1];
P = polytope(A,b);
res(end+1,1) = ~isFullDim(P) && ~isempty(P.fullDim.val) && ~P.fullDim.val;

% 3D, non-degenerate, unbounded (with subspace computation)
% ...unit square in x1-x2, unbounded in x3
A = [1 0 0;-1 0 0;0 1 0;0 -1 0]; b = [1;1;1;1];
P = polytope(A,b);
[res_, X] = isFullDim(P);
res(end+1,1) = res_ && rank(X) == 3 && all(size(X) == [3,3]);


% 4D, degenerate, unbounded (with subspace computation)
% ...unit square in x1-x2, unbounded in x3, x4 = 0
A = [1 0 0 0;-1 0 0 0;0 1 0 0;0 -1 0 0]; b = [1;1;1;1];
Ae = [0 0 0 1]; be = 0;
P = polytope(A,b,Ae,be);
[res_, X] = isFullDim(P);
res(end+1,1) = ~res_ && rank(X) == 3 && all(size(X) == [4,3]);
% same polytope, but represented via only inequality constraints
A = [1 0 0 0;-1 0 0 0;0 1 0 0;0 -1 0 0; 0 0 0 1; 0 0 0 -1];
b = [1;1;1;1;0;0];
P = polytope(A,b);
[res_, X] = isFullDim(P);
res(end+1,1) = ~res_ && rank(X) == 3 && size(X,1) == 4 && size(X,2) == 3;

% 4D, degenerate, bounded (with subspace computation)
% ...unit square in x1-x2, x3 = 0
A = [1 0 0;-1 0 0;0 1 0;0 -1 0]; b = [1;1;1;1]; Ae = [0 0 1]; be = 0;
P = polytope(A,b,Ae,be);
[res_, X] = isFullDim(P);
res(end+1,1) = ~res_ && rank(X) == 2 && all(size(X) == [3,2]);
% same polytope, but represented via only inequality constraints
A = [1 0 0;-1 0 0;0 1 0;0 -1 0;0 0 1;0 0 -1]; b = [1;1;1;1;0;0];
P = polytope(A,b);
[res_, X] = isFullDim(P);
res(end+1,1) = ~res_ && rank(X) == 2 && all(size(X) == [3,2]);


% sequence of functions
% 2D, full-dimensional
A = [1 1; -2 1; -4 -2; 2 -3]; b = ones(4,1);
P = polytope(A,b);
isFullDim(P);

% intersection of 2 non-degenarate polytopes could be both degenarate or not -> empty property
A = [1 0; -1 0; 0 1; 0 -1]; b = [1;1;1;1];
P2 = polytope(A,b);
res(end+1,1) = isFullDim(P2);
P = P & P2;
res(end+1,1) = isempty(P.fullDim.val);

% intersection with degenerate polytope would always be degenarate
A = [1 0; -1 0; 0 1; 0 -1]; b = [1;0;0;0];
P3 = polytope(A,b);
res(end+1,1) = ~isFullDim(P3);
P = P & P3;
res(end+1,1) = ~isempty(P.fullDim.val) && ~P.fullDim.val;

% no change
P = normalizeConstraints(P);
res(end+1,1) = ~isempty(P.fullDim.val) && ~P.fullDim.val;

% Minkowski sum with fully dimensional polytope would be fully dimensional
A = [0.9691, -0.2466; -0.2466, -0.9691; -0.5109, 0.8597; 0.6113, -0.7914; -0.4098, 0.9122];
b = [0.9185; 1.0376; 1.0311; 0.9590; 1.0213];
P4 = polytope(A,b);
res(end+1,1) = isFullDim(P4);
P = P + P4;
res(end+1,1) = ~isempty(P.fullDim.val) && P.fullDim.val;

% Projecting to higher dimension changes fullDim property
P_ = projectHighDim(P,4,[3,4]);
res(end+1,1) = ~isempty(P_.fullDim.val) && ~P_.fullDim.val;

% Lifting to higher dimension does not change fullDim property
P = lift(P,4,[3,4]);
res(end+1,1) = ~isempty(P.fullDim.val) && P.fullDim.val;

% If one polytope is non-degenerate, then the Minkowski diff is non-degenerate
A = [
    0.3774, -0.6262, -0.6131, -0.2993 ; ...
    -0.3539, -0.3854, 0.5092, -0.6833 ; ...
    0.3673, -0.4032, 0.0985, 0.8323 ; ...
    0.5925, 0.6134, 0.3778, 0.3605 ; ...
    -0.6365, 0.3594, -0.6019, -0.3216 ; ...
    -0.6572, -0.3818, -0.6065, 0.2336 ; ...
    0.3430, 0.8561, 0.2134, -0.3225 ; ...
    0.1438, -0.7332, 0.5918, -0.3025 ; ...
    -0.6270, -0.4834, -0.4739, -0.3856 ; ...
    -0.1397, 0.4827, 0.3141, 0.8055 ; ...
    0.5720, 0.5746, 0.0499, -0.5832 ; ...
    -0.5720, -0.5746, -0.0499, 0.5832 ; ...
 ];
b = [1.6556; 0.2467; 2.7103; 1.1828; -0.1755; 1.2791; -0.0778; 1.4433; 0.6343; 1.3312; -1.9076; 1.9076];
P5 = polytope(A,b);
P = minkDiff(P,P5);
res(end+1,1) = ~isempty(P.fullDim.val) && P.fullDim.val;

% linear map with invertible matrix keeps properties
A = [1   3  -2   4;
     0   1   2  -1;
    -1   2   0   3;
     1  -1   3   2;];
P = A*P;
res(end+1,1) = ~isempty(P.fullDim.val) && P.fullDim.val;

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
