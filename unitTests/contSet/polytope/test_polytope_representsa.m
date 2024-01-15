function res = test_polytope_representsa
% test_polytope_representsa - unit test function of representsa
%
% Syntax:
%    res = test_polytope_representsa
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

% Authors:       Mark Wetzlinger, Viktor Kotsev
% Written:       17-March-2022
% Last update:   16-December-2023 (MW, more comparisons)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true(0);

% --- origin --------------------------------------------------------------

% fully empty case
P = polytope.empty(2);
res(end+1,1) = ~representsa(P,'origin');

% 2D, fully empty
A = zeros(0,2); b = zeros(0,0);
P = polytope(A,b);
res(end+1,1) = ~representsa(P,'origin');

% 2D, only origin
A = [1 0; 0 1; -1 0; 0 -1]; b = zeros(4,1);
P = polytope(A,b);
res(end+1,1) = representsa(P,'origin');

% 2D, shifted center
P = P + [0.01; 0];
res(end+1,1) = ~representsa(P,'origin');
% ...add tolerance
tol = 0.02;
res(end+1,1) = representsa(P,'origin',tol);


% --- interval ------------------------------------------------------------

% 2D, fully empty
A = zeros(0,2); b = zeros(0,0);
P = polytope(A,b);
[res(end+1,1),I] = representsa(P,'interval');
I_true = interval(-Inf(2,1),Inf(2,1));
res(end+1,1) = isequal(I,I_true);

% 2D, not an interval
A = [2 1; 1 0; 0 1]; b = [1; 1; 2];
P = polytope(A,b);
res(end+1,1) = ~representsa(P,'interval');

% 3D, unit box
n = 3; A = [eye(n); -eye(n)]; b = ones(2*n,1);
P = polytope(A,b);
[res(end+1,1),I] = representsa(P,'interval');
I_true = interval(-ones(n,1),ones(n,1));
res(end+1,1) = isequal(I,I_true);

% 3D, degenerate interval
% A = [1 0 0; -1 0 0]; b = [1;1]; Ae = [0 1 0; 0 0 1]; be = [5; 7];
% P = polytope(A,b,Ae,be);
% [res(end+1,1),I] = representsa(P,'interval');
% res(end+1,1) = isequal(I,interval([-1;5;7],[1;5;7]));

% 3D, bounded
n = 3;
A = [2; 3; 2.75; 4; 4.5; 1.5] .* [eye(n); -eye(n)];
b = [0.5; 4; 0.25; 2; 2.25; 3] .* [ones(2*n,1)];
P = polytope(A,b);
[res(end+1,1),I] = representsa(P,'interval');
lb_true = [-0.5;-0.5;-2]; ub_true = [0.25;4/3;1/11];
I_true = interval(lb_true,ub_true);
res(end+1,1) = isequal(I,I_true);

% 3D, unbounded polytope (still an interval)
n = 3;
A = [2*eye(n); -eye(n)]; A([2,4],:) = [];
b = [ones(2*n,1)]; b([2,4]) = [];
P = polytope(A,b);
[res(end+1,1),I] = representsa(P,'interval');
lb_true = [-Inf;-1;-1]; ub_true = [0.5;Inf;0.5];
I_true = interval(lb_true,ub_true);
res(end+1,1) = isequal(I,I_true);


% --- emptySet ------------------------------------------------------------

% empty constructor
P = polytope.empty(2);
res(end+1,1) = representsa(P,'emptySet');

% 1D, fully empty
A = zeros(0,1); b = zeros(0,0);
P = polytope(A,b);
res(end+1,1) = ~representsa(P,'emptySet');

% 1D, empty
A = [1; -1]; b = [1; -3];
P = polytope(A,b);
res(end+1,1) = representsa(P,'emptySet');

% 1D, empty
A = [1; -1; 1; 1; -1; -1]; b = [1; -3; 1; 4; 2; 1];
P = polytope(A,b);
res(end+1,1) = representsa(P,'emptySet');


% 2D, polytope encloses origin
A = [2 1; -2 3; -2 -2; 4 1]; b = ones(4,1);
P = polytope(A,b);
res(end+1,1) = ~representsa(P,'emptySet');

% 2D, empty, only inequalities
A = [-1 -1; -1 1; 1 0]; b = [-2; -2; -1];
P = polytope(A,b);
res(end+1,1) = representsa(P,'emptySet');

% 2D, empty only equalities: x1 == 1, x2 == 1, x1+x2 == 1
Ae = [1 0; 0 1; 1 1]; be = [1;1;1];
P = polytope([],[],Ae,be);
res(end+1,1) = representsa(P,'emptySet');

% 2D, did not work using MPT toolbox (both should be non-empty!)
A = [0.000, 0.707, -0.707; ...
    -0.707, 0.000, 0.707; ...
    0.707, -0.707, 0.000; ...
    0.000, 0.000, 1.000; ...
    0.000, -1.000, 0.000; ...
    1.000, 0.000, 0.000; ...
    0.000, -0.707, 0.707; ...
    0.707, 0.000, -0.707; ...
    -0.707, 0.707, 0.000; ...
    0.000, 0.000, -1.000; ...
    0.000, 1.000, 0.000; ...
    -1.000, 0.000, 0.000];
b1 = [0.9306; 0.9306; 0.7693; 1.7720; 1.5440; 1.5440; 0.9306; 0.9306; ...
      0.7693; 1.7720; 1.5440; 1.5440];
b2 = [0.9284; 0.9284; 0.7665; 1.7710; 1.5420; 1.5420; 0.9284; 0.9284; ...
      0.7665; 1.7710; 1.5420; 1.5420];
P1 = polytope(A,b1);
P2 = polytope(A,b2);
res(end+1,1) = ~representsa(P1,'emptySet');
res(end+1,1) = ~representsa(P2,'emptySet');

% 2D, unbounded, degenerate (line)
A = [0 1;0 -1]; b = [1;-1];
P = polytope(A,b);
res(end+1,1) = ~representsa(P,'emptySet');

% 2D, V-representation
V = [2 0; -2 0; 0 2; 0 -2]';
P = polytope(V);
res(end+1,1) = ~representsa(P,'emptySet');


% --- conHyperplane -------------------------------------------------------

% 2D, degenerate
A = [1 1;1 0;-1 -1]; b = [1;2;-1];
P = polytope(A,b);
[res(end+1,1),hyp] = representsa(P,'conHyperplane');
% init equivalent conHyperplane
hyp_ = conHyperplane(halfspace([0.5 0.5],0.5),[1 0],2);
res(end+1,1) = isequal(hyp,hyp_);

% 2D, degenerate
A = [-1 0; 0 1]; b = [0;0]; Ae = [1 0]; be = 0;
P = polytope(A,b,Ae,be);
[res(end+1,1),hyp] = representsa(P,'conHyperplane');
% init equivalent conHyperplane
hyp_ = conHyperplane(halfspace([1 0],0),[-1 0; 0 1],[0;0]);
res(end+1,1) = isequal(hyp,hyp_);


% --- point ---------------------------------------------------------------

% 2D, fully empty
A = zeros(0,2); b = zeros(0,0);
P = polytope(A,b);
res(end+1,1) = ~representsa(P,'point');

% 5D, only equality constraints
n = 5; [Ae,~,~] = svd(randn(n)); be = ones(n,1);
P = polytope([],[],Ae',be);
% compute box and check if its a point
B = box(P);
res(end+1,1) = representsa(B,'point');


% --- fullspace -----------------------------------------------------------

% 2D, fully empty
A = zeros(0,2); b = zeros(0,0);
P = polytope(A,b);
[res(end+1,1),fs] = representsa(P,'fullspace');
res(end+1,1) = isequal(fs,fullspace(2));

% 2D, Inf
P = polytope.Inf(2);
[res(end+1,1),fs] = representsa(P,'fullspace');
res(end+1,1) = isequal(fs,fullspace(2));

% 3D
A = [0,0,0]; b = 4;
P = polytope(A,b);
[res(end+1,1),fs] = representsa(P,'fullspace');
res(end+1,1) = isequal(fs,fullspace(3));


% --- conZonotope ---------------------------------------------------------

% 2D, fully empty
A = zeros(0,2); b = zeros(0,0);
P = polytope(A,b);
res(end+1,1) = ~representsa(P,'conZonotope');

% 2D, bounded
A = [1 0; -1 1; -1 -1]; b = [1;1;1];
P = polytope(A,b);
[res(end+1,1),cZ] = representsa(P,'conZonotope');
res(end+1,1) = isequal(cZ,P);

% 3D, degenerate
A = [1 0 0; -1 1 0; -1 -1 0]; b = [1;1;1]; Ae = [0 0 1]; be = 0;
P = polytope(A,b,Ae,be);
[res(end+1,1),cZ] = representsa(P,'conZonotope');
res(end+1,1) = isequal(cZ,P);

% 2D, unbounded
A = [1 0; 0 1]; b = [1; 2];
P = polytope(A,b);
res(end+1,1) = ~representsa(P,'conZonotope');


% --- ellipsoid -----------------------------------------------------------

% 1D, bounded
A = [1;-1]; b = [5;-2];
P = polytope(A,b);
[res(end+1,1),E] = representsa(P,'ellipsoid');
res(end+1,1) = compareMatrices(vertices(P),vertices(E),tol);

% 1D, unbounded
A = 1; b = 1;
P = polytope(A,b);
res(end+1,1) = ~representsa(P,'ellipsoid');

% 2D, bounded
A = [1 0; -1 1; -1 -1]; b = [1;1;1];
P = polytope(A,b);
res(end+1,1) = ~representsa(P,'ellipsoid');

% 2D, degenerate (point)
Ae = [1 0; 0 1]; be = [2; 3];
P = polytope([],[],Ae,be);
[res(end+1,1),E] = representsa(P,'ellipsoid');
res(end+1,1) = all(withinTol(E.q,[2;3]));

% 2D, unbounded
A = [1 0]; b = 2;
P = polytope(A,b);
res(end+1,1) = ~representsa(P,'ellipsoid');

% 2D, empty
% P = polytope.empty(2);
% res(end+1,1) = ~representsa(P,'ellipsoid');

% 2D, degenerate (bounded line)
% A = [1 0; -1 1; -1 1]; b = [1;1;1]; Ae = [1 0]; be = 0;
% P = polytope(A,b,Ae,be);
% res(end+1,1) = representsa(P,'ellipsoid');


% --- polyZonotope --------------------------------------------------------

% 1D, bounded
A = [1;-1]; b = [5;-2];
P = polytope(A,b);
[res(end+1,1),pZ] = representsa(P,'polyZonotope');
res(end+1,1) = withinTol(supportFunc_(pZ,1,'upper','globOpt',6,1e-6),5,1e-6) ...
    && withinTol(supportFunc_(pZ,1,'lower','globOpt',6,1e-6),2,1e-6);

% 1D, unbounded
A = 1; b = 1;
P = polytope(A,b);
res(end+1,1) = ~representsa(P,'polyZonotope');

% 2D, bounded
A = [1 0; -1 1; -1 -1]; b = [1;1;1];
P = polytope(A,b);
res(end+1,1) = representsa(P,'polyZonotope');

% 2D, degenerate (point)
Ae = [1 0; 0 1]; be = [2; 3];
P = polytope([],[],Ae,be);
[res(end+1,1),pZ] = representsa(P,'polyZonotope');
res(end+1,1) = all(withinTol(pZ.c,[2;3]));

% 2D, unbounded
A = [1 0]; b = 2;
P = polytope(A,b);
res(end+1,1) = ~representsa(P,'polyZonotope');


% --- zonoBundle ----------------------------------------------------------

% 1D, unbounded
A = 1; b = 1;
P = polytope(A,b);
res(end+1,1) = ~representsa(P,'zonoBundle');

% 2D, unbounded
A = [4 1; 2 3; 0 2; -1 4; -3 1]; b = [1;1;1;1;1];
P = polytope(A,b);
res(end+1,1) = ~representsa(P,'zonoBundle');

% 2D, bounded
A = [1 0; -1 1; -1 -1]; b = [1;1;1];
P = polytope(A,b);
[res(end+1,1),zB] = representsa(P,'zonoBundle');
res(end+1,1) = isequal(P,zB,1e-10);


% --- halfspace -----------------------------------------------------------

% 1D, unbounded
A = -1; b = 1;
P = polytope(A,b);
[res(end+1,1),hs] = representsa(P,'halfspace');
res(end+1,1) = withinTol(hs.c,-1) && withinTol(hs.d,1);

% 1D, bounded
A = [1;-1]; b = [5;-2];
P = polytope(A,b);
res(end+1,1) = ~representsa(P,'halfspace');

% 2D, unbounded with redundant halfspaces
A = [1 0; 1 0; 1 0]; b = [4;3;2];
P = polytope(A,b);
[res(end+1,1),hs] = representsa(P,'halfspace');
res(end+1,1) = all(withinTol(hs.c,[1;0])) && withinTol(hs.d,2);


% --- conPolyZono ---------------------------------------------------------

% 1D, unbounded
A = -1; b = 1;
P = polytope(A,b);
res(end+1,1) = ~representsa(P,'conPolyZono');

% 2D, degenerate, bounded
A = [0 1; 0 -1]; b = [1;1]; Ae = [1 0]; be = 2;
P = polytope(A,b,Ae,be);
[res(end+1,1),cPZ] = representsa(P,'conPolyZono');
res(end+1,1) = withinTol(supportFunc(cPZ,[1;0],'upper'),2,1e-10) ...
    && withinTol(supportFunc(cPZ,[0;1],'upper'),1,1e-10) ...
    && withinTol(supportFunc(cPZ,[0;1],'lower'),-1,1e-10);


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
