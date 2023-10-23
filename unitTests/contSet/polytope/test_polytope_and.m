function res = test_polytope_and
% test_polytope_and - unit test function of intersection between a polytope
%    and another set or vector
%
% Syntax:
%    res = test_polytope_and
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

% Authors:       Viktor Kotsev
% Written:       09-May-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% TODOs
% intersection with fully-empty polytope
% equality constraints
% unbounded?


% 1D (convertible to intervals) -------------------------------------------

% empty intersection
% P1 ... x <= 2, x >= 1
P1 = polytope([1; -1],[2;-1]);
% P2 ... x <= 5, x >= 4
P2 = polytope([1;-1],[5; -4]);
% compute intersection and check emptiness
P_ = P1 & P2;
res = representsa(P_,'emptySet');


% intersection at only one point
% P1 ... x <= 2, x >= 1
P1 = polytope([1; -1],[2;-1]);
% P2 ... x <= 3, x >= 2
P2 = polytope([1;-1],[3;-2]);
% check boundedness for property test (should be true)
isBounded(P1);
isBounded(P2);
% compute intersection, compare to true result (by hand) and check
% boundedness
P_ = P1 & P2;
Ptrue = polytope([1;-1], [2;-2]);
res(end+1,1) = (P_ == Ptrue) && ~isempty(P_.bounded.val) && P_.bounded.val;


% intersection is a full-dimensional polytope
% P1 ... x <= 3, x >= 1
P1 = polytope([1;-1], [3;-1]);
% P2 ... x <= 4, x >= 2
P2 = polytope([1;-1], [4;-2]);
% compute intersection and compare to true result
P_ = P1 & P2;
Ptrue = polytope([1;-1],[3;-2]);
res(end+1,1) = P_ == Ptrue;

% emptiness property test
% P1 ... x <= 2, x >= 3 (empty)
P1 = polytope([1; -1],[2;-3]);
% P2 ... x <= 3, x >= 2
P2 = polytope([1;-1],[3;-2]);
% determine emptiness of P1
representsa(P1,"emptySet");
% compute intersection and check emptiness
P_ = P1 & P2;
res(end+1,1) = ~isempty(P_.emptySet.val) && P_.emptySet.val;


% 2D ----------------------------------------------------------------------

% bounded & unbounded -> bounded
A = [-1 1; 0 1; 1 0; -1 -2];
b = [0; 1; 2; 0];
P1 = polytope(A,b);
% determine boundedness
isBounded(P1);
% unbounded polytope
P2 = polytope([-1 -1]/sqrt(2),-sqrt(2));
% compute intersection, compare to true result (by hand) and check
% boundedness
P_ = P1 & P2;
Ptrue = polytope([-1/sqrt(2) -1/sqrt(2); 0 1; 1 0],[-sqrt(2); 1; 2]);
res(end+1,1) = P_ == Ptrue && ~isempty(P_.bounded.val) && P_.bounded.val;

% visualization for debugging
% figure; hold on; axis([-3,3,-3,3]);
% plot(halfspace([-1 -1]/sqrt(2),-sqrt(2)));
% plot(P1,[1,2],'k');
% plot(Ptrue,[1,2],'g');
% plot(P_,[1,2],'r--');
% close;

% bounded & unbounded -> bounded
A = [-1 1; 0 1; 1 0; -1 -2];
b = [0; 1; 2; 0];
P1 = polytope(A,b);
% determine boundedness
isBounded(P1);
% unbounded polytope
P2 = polytope([-1 2],0);
% compute intersection and compare to true result
P_ = P1 & P2;
Ptrue = polytope([-1 2; -1 -2; 1 0],[0; 0; 2]);
res(end+1,1) = P_ == Ptrue;

% visualization for debugging
% figure; hold on; axis([-3,3,-3,3]);
% plot(halfspace([-1 2],0));
% plot(P1,[1,2],'k');
% plot(Ptrue,[1,2],'g');
% plot(P_,[1,2],'r--');
% close;


% bounded & bounded -> bounded
A = [-1 1; 0 1; 1 0; -1 -2];
b = [0; 1; 2; 0];
P1 = polytope(A,b);
% determine boundedness
isBounded(P1);
% another bounded polytope
P2 = polytope([-1/sqrt(2) -1/sqrt(2); -1 2],[-sqrt(2); 0]);
% compute intersection, compare to true result
P_ = P1 & P2;
Ptrue = polytope([-1/sqrt(2) -1/sqrt(2); -1 2; 1 0],[-sqrt(2); 0; 2]);
% intersection
res(end+1,1) = P_ == Ptrue;

% visualization
% figure; hold on; axis([-3,3,-3,3]);
% plot(halfspace([-1/sqrt(2) -1/sqrt(2)],-sqrt(2)));
% plot(halfspace([-1 2],0),[1,2],'m');
% plot(P1,[1,2],'k');
% plot(Ptrue,[1,2],'g');
% plot(P_,[1,2],'r--');
% close;

% unbounded & unbounded -> bounded
P1 = polytope([1 1; 1 -1],[2; 2]);
P2 = polytope([-1 1;-1 -1],[2; 2]);
P_ = P1 & P2;
Ptrue = polytope([1 1; 1 -1;-1 1;-1 -1], [2; 2; 2; 2]);
res(end+1,1) = P_ == Ptrue;


% empty & bounded -> empty
A1 = [-1 -1; -1 1; 1 0];
b1 = [-2; -2; -1];
P1 = polytope(A1,b1);
A2 = [-1 1; 0 1; 1 0; -1 -2];
b2 = [0; 1; 2; 0];
P2 = polytope(A2,b2);
% determine emptiness of P1
representsa(P1,"emptySet");
% compute intersection and check emptiness value
P_ = P1 & P2;
res(end+1,1) = ~isempty(P_.emptySet.val) && P_.emptySet.val;


% nD ----------------------------------------------------------------------

% 3D: degenerate & non-degenerate
% square in 3D
A1 = [0 0 1; 0 0 -1; 1 1 0; -1 1 0; 1 -1 0; -1 -1 0];
b1 = [1;-1;2;2;2;2];
P1 = polytope(A1,b1);
% obtain information about degeneracy
isFullDim(P1);
% init non-degenerate polytope
A2 = [0 0 1; 0 0 -1; 1 1 0; -1 1 0; 1 -1 0; -1 -1 0];
b2 = [1;-1;3;-1;3;-1];
P2 = polytope(A2,b2);
% compute intersection
P_ = P1 & P2;
% true solution (by hand)
Ptrue = polytope([0 0 1; 0 0 -1; 1 1 0; -1 1 0; 1 -1 0; -1 -1 0],[1; -1; 2; -1; 2; -1]);
% compare, also P_ should be degenerate since P1 is degenerate
res(end+1,1) = P_ == Ptrue && ~isempty(P_.fullDim.val) && ~P_.fullDim.val;


% combine result
res = all(res);

% ------------------------------ END OF CODE ------------------------------
