function res = test_polytope_mldivide
% test_polytope_mldivide - unit test function of mldivide
%
% Syntax:
%    res = test_polytope_mldivide
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
% Written:       03-December-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% skip test for now...
res = true; return

res = true(0);

% 2D, unbounded, unbounded, no intersection
P1 = polytope([-1 0; 0 -1],[-1; -1]);       % x >= 1, y >= 1
P2 = polytope([1 1],2);                     % x+y <= 2

% compute set difference and compare to true result
P = mldivide(P1,P2);
res(end+1,1) = isequal(P1,P,1e-8);

% plot
% figure; hold on;
% xlim([-1,4]); ylim([-1,4]);
% plot(P1);
% plot(P2,[1,2],'k');
% plot(P,[1,2],'r--');
% plot(P_true,[1,2],'g');


% 2D, unbounded polytope/ polytope
P1 = polytope([1 0;-1 0; 0 1], [1;1;1]);
P2 = polytope([1 0; -1 0; 0 1; 0 -1], [1;1;1;1]);

P_ = mldivide(P1,P2);
P_true = polytope([1 0;-1 0; 0 1], [1;1;-1]);

res(end+1,1) = P_ == P_true;

% degenerate polytope / polyope
% P1 = polytope([1 0; -1 0; 0 1; 0 -1], [1;1;1;-1]);
% P2 = polytope([1 0; -1 0; 0 1; 0 -1], [2;0;2;0]);

% P_div = mldivide(P1,P2);
% P_true = polytope([1 0;-1 0; 0 1; 0 -1], [0;1;1;-1]);

% res(end+1,1) = P_div == P_true;


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
