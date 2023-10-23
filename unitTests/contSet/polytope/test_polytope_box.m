function res = test_polytope_box
% test_polytope_box - unit test function of box
%
% Syntax:
%    res = test_polytope_box
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
% Written:       16-May-2022
% Last update:   27-July-2023 (MW, more tests)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true(0);

% 1D ----------------------------------------------------------------------
% all box conversions are exact since 1D polytope are intervals

% bounded
P = polytope([1;-1],[1;1]);
B = box(P);
res(end+1,1) = B == P;

% unbounded
P = polytope(3,6);
B = box(P);
res(end+1,1) = B == P;

% single point (inequalities)
P = polytope([2;-1],[1;-0.5]);
B = box(P);
res(end+1,1) = B == P;

% single point (equality)
P = polytope([],[],1,3);
B = box(P);
res(end+1,1) = B == P;

% empty (inequalities)
P = polytope([2;-1],[1;-4]);
B = box(P);
res(end+1,1) = isemptyobject(B);

% empty (equalities)
P = polytope([],[],[2;1],[5;2]);
B = box(P);
res(end+1,1) = isemptyobject(B);


% 2D ----------------------------------------------------------------------

% bounded 
A = [1 1; 1 -1;-1 1;-1 -1];
b = [2; 2; 2; 2];
P = polytope(A,b);
% compute box and compare to true result
B = box(P);
Ptrue = polytope([1 0; 0 1;-1 0; 0 -1],[2; 2; 2; 2]);
res(end+1,1) = B == Ptrue;

% bounded and degenerate
A = [1 0; -2 0]; b = [2; 5];
Ae = [-1 1]; be = 2;
P = polytope(A,b,Ae,be);
% compute box and compare to true result
B = box(P);
Ptrue = polytope([1 0; 0 1;-1 0; 0 -1],[2; 4; 2.5; 0.5]);
res(end+1,1) = B == Ptrue;

% unbounded
P = polytope([1 1; 1 -1],[2; 2]);
B = box(P);
% compute box and compare to true result
Ptrue = polytope([1 0],2);
res(end+1,1) = Ptrue == B;

% unbounded and degenerate
P = polytope([],[],[1 1],2);
B = box(P);
res(end+1,1) = isemptyobject(B);

% unbounded and degenerate (axis-aligned)
P = polytope([],[],[0 3],2);
B = box(P);
% compute box and compare to true result
Ptrue = polytope([],[],[0 1],2/3);
res(end+1,1) = Ptrue == B;

% bounded by vertex instantiation
V = [1 0; -1 0; 0 1; 0 -1]';
P = polytope(V);
% compute box and compare to true result
B = box(P);
Ptrue = polytope([1 0;-1 0;0 1;0 -1],[1;1;1;1]);
res(end+1,1) = Ptrue == B;


% nD ----------------------------------------------------------------------

% 3D, degenerate
A = [0 0 1; 0 0 -1; 1 1 0; -1 1 0; 1 -1 0; -1 -1 0];
b = [1;-1;2;2;2;2];
P = polytope(A,b);
% compute box and compute to true result
B = box(P);
Ptrue = polytope([0 0 1; 0 0 -1; 1 0 0; -1 0 0; 0 1 0; 0 -1 0], [1;-1;2;2;2;2]);
res(end+1,1) = Ptrue == B;

% 5D, single point
[Ae,~,~] = svd(randn(5));
be = ones(5,1);
P = polytope([],[],Ae',be);
% compute box and check if its a point
B = box(P);
res(end+1,1) = representsa(B,'point');


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
