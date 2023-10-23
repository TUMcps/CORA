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
% See also: -

% Authors:       Mark Wetzlinger, Viktor Kotsev
% Written:       17-March-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true(0);

% 1. compare to origin
% empty case
res(end+1,1) = ~representsa(polytope(),'origin');

% only origin
P = polytope([1 0; 0 1; -1 0; 0 -1],zeros(4,1));
res(end+1,1) = representsa(P,'origin');

% shifted center
P = P + [0.01; 0];
res(end+1,1) = ~representsa(P,'origin');
% add tolerance
tol = 0.02;
res(end+1,1) = representsa(P,'origin',tol);


% 2. compare to interval
% unit box
P = polytope([eye(3); -eye(3)],ones(6,1));
[res(end+1,1),I] = representsa(P,'interval');
res(end+1,1) = isequal(I,interval(-ones(3,1),ones(3,1)));

% degenerate interval
% P = polytope([1 0 0; -1 0 0],[1;1],[0 1 0; 0 0 1],[5; 7]);
% [res(end+1,1),I] = representsa(P,'interval');
% res(end+1,1) = isequal(I,interval([-1;5;7],[1;5;7]));

% init unit box as a polytope
n = 3;
A = [eye(n); -eye(n)]; b = [ones(2*n,1)];
P = polytope(A,b);

% check if polytope is an interval and compute corresponding interval
[res(end+1,1),I] = representsa(P,'interval');
res(end+1,1) = isequal(I,interval(-ones(n,1),ones(n,1)));


% use different values in A and b
n = 3;
A = [2; 3; 2.75; 4; 4.5; 1.5] .* [eye(n); -eye(n)];
b = [0.5; 4; 0.25; 2; 2.25; 3] .* [ones(2*n,1)];
P = polytope(A,b);
[res(end+1,1),I] = representsa(P,'interval');
I_true = interval([-0.5;-0.5;-2],[0.25;4/3;1/11]);
res(end+1,1) = isequal(I,I_true);


% unbounded polytope (still an interval)
n = 3;
A = [2*eye(n); -eye(n)];
A([2,4],:) = [];
b = [ones(2*n,1)];
b([2,4]) = [];
P = polytope(A,b);

[res(end+1,1),I] = representsa(P,'interval');
I_true = interval([-Inf;-1;-1],[0.5;Inf;0.5]);
res(end+1,1) = isequal(I,I_true);


% not an interval
A = [2 1; 1 0; 0 1]; b = [1; 1; 2];
P = polytope(A,b);
res(end+1,1) = ~representsa(P,'interval');


% 3. compare to empty set
% fully empty case
res(end+1,1) = representsa(polytope(),'emptySet');

% 1D cases
% init polytope where two halfspaces cannot be fulfilled at the same time:
%   x <= 1, x >= 3
A = [1; -1]; b = [1; -3];
P = polytope(A,b);
res(end+1,1) = representsa(P,'emptySet');

% add other halfspaces (which obviously do not change the result)
A = [1; -1; 1; 1; -1; -1]; b = [1; -3; 1; 4; 2; 1];
P = polytope(A,b);
res(end+1,1) = representsa(P,'emptySet');


% 2D cases
% init polytope enclosing origin
A = [2 1; -2 3; -2 -2; 4 1];
b = ones(4,1);
P = polytope(A,b);
res(end+1,1) = ~representsa(P,'emptySet');

% 2D, only inequalities, empty
%   x1 + x2 >= 2, x1 - x2 >= 2, x1 <= -1
A = [-1 -1; -1 1; 1 0]; b = [-2; -2; -1];
P = polytope(A,b);
res(end+1,1) = representsa(P,'emptySet');
% visualization
% figure; hold on; axis([-4,4,-4,4]);
% plot(halfspace(A(1,:),b(1)));
% plot(halfspace(A(2,:),b(2)),[1,2],'FaceColor',colorblind('r'));
% plot(halfspace(A(3,:),b(3)),[1,2],'FaceColor',colorblind('y'));


% 2D, only equalities, empty
%    x1 == 1, x2 == 1, x1+x2 == 1
P = polytope([],[],[1 0; 0 1; 1 1],[1;1;1]);
res(end+1,1) = representsa(P,'emptySet');


% test below did not work using MPT toolbox
% define values for which mpt toolbox outputs wrong results (both should be
% non-empty!)
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

b1 = [0.9306
     0.9306
     0.7693
     1.7720
     1.5440
     1.5440
     0.9306
     0.9306
     0.7693
     1.7720
     1.5440
     1.5440];
b2 = [0.9284
     0.9284
     0.7665
     1.7710
     1.5420
     1.5420
     0.9284
     0.9284
     0.7665
     1.7710
     1.5420
     1.5420];

% instantiate polytopes
P1 = polytope(A,b1);
P2 = polytope(A,b2);

% check results
res(end+1,1) = ~representsa(P1,'emptySet');
res(end+1,1) = ~representsa(P2,'emptySet');

% unbounded degenerate 2D line
A = [0 1;0 -1]; b = [1;-1];
P = polytope(A,b);
res(end+1,1) = ~representsa(P,'emptySet');

% V-representation test
V = [2 0; -2 0; 0 2; 0 -2]';
P = polytope(V);
res(end+1,1) = ~representsa(P,'emptySet');


% 4. compare to constrained hyperplane
A = [1 1;1 0;-1 -1]; b = [1;2;-1];
P = polytope(A,b);
[res(end+1,1),hyp] = representsa(P,'conHyperplane');
% init equivalent conHyperplane
hyp_ = conHyperplane(halfspace([0.5 0.5],0.5),[1 0],2);
res(end+1,1) = isequal(hyp,hyp_);

A = [-1 0; 0 1]; b = [0;0];
Ae = [1 0]; be = 0;
P = polytope(A,b,Ae,be);
[res(end+1,1),hyp] = representsa(P,'conHyperplane');
% init equivalent conHyperplane
hyp_ = conHyperplane(halfspace([1 0],0),[-1 0; 0 1],[0;0]);
res(end+1,1) = isequal(hyp,hyp_);

% 5. compare to point
[Ae,~,~] = svd(randn(5));
be = ones(5,1);
P = polytope([],[],Ae',be);
% compute box and check if its a point
B = box(P);
res(end+1,1) = representsa(B,'point');

% 6. compare to fullspace
P = polytope([0,0,0],4);
[res(end+1,1),fs] = representsa(P,'fullspace');
res(end+1,1) = isequal(fs,fullspace(3));


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
