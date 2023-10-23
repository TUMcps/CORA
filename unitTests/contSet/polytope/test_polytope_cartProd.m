function res = test_polytope_cartProd
% test_polytope_cartProd - unit test function of Cartesian product
%
% Syntax:
%    res = test_polytope_cartProd
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
% Last update:   27-July-2023
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true(0);

% 1D and 2D (both only inequalities)
P1 = polytope(1,0);
P2 = polytope([1 1; -1 -1],[1;1]);
P_ = cartProd(P1,P2);
A_true = [1, 0, 0; 0, 1, 1; 0,-1,-1;];
b_true = [0;1;1];
P_true = polytope(A_true,b_true);
% check for equality
res(end+1,1) = isequal(P_,P_true,1e-10);


% 2D and 2D (both only inequalities)
P1 = polytope([-1 -1; 1 0;-1 0; 0 1; 0 -1],[2; 3; 2; 3; 2]);
P2 = polytope([1;-1],[3; 2]);
P_ = cartProd(P1,P2);
A_true = [-1, -1, 0; ...
           1,  0, 0; ...
          -1,  0, 0; ...
           0,  1, 0; ...
           0, -1, 0; ...
           0,  0, 1; ...
           0,  0,-1;];
b_true = [2;3;2;3;2;3;2];
P_true = polytope(A_true,b_true);
% check for equality and emptiness
res(end+1,1) = isequal(P_,P_true,1e-10);


% 1D and 3D (both only equalities)
P1 = polytope([],[],5,1);
P2 = polytope([],[],[1 0 1; -1 0 0; 0 1 1],[2;4;-3]);
P_ = cartProd(P1,P2);
A_true = [1 0 0 0; 0 1 0 1; 0 -1 0 0; 0 0 1 1];
b_true = [0.2;2;4;-3];
P_true = polytope(zeros(0,4),[],A_true,b_true);
% check for equality
res(end+1,1) = isequal(P_,P_true,1e-10);


% 2D (only equalities) and 3D (only inequalities)
P1 = polytope(zeros(0,2),[],[1 0; 0 1],[1;1]);
P2 = polytope([1 0 0; 0 1 0; 0 0 1; -1 -1 -1],ones(4,1));
P_ = cartProd(P1,P2);
A_true = [0 0 1 0 0; 0 0 0 1 0; 0 0 0 0 1; 0 0 -1 -1 -1];
b_true = ones(4,1);
Ae_true = [1 0 0 0 0; 0 1 0 0 0];
be_true = [1;1];
P_true = polytope(A_true,b_true,Ae_true,be_true);
% check for equality
res(end+1,1) = isequal(P_,P_true,1e-10);


% 1D and 4D (both with equalities and inequalities)
P1 = polytope([1;2],[4;20],1,0); % with redundancies
P2 = polytope([2 1 0 0; -1 1 0 0; 0 -4 0 0],[3;2;5],[0 0 1 1],1);
P_ = cartProd(P1,P2);
A_true = [0 2 1 0 0; 0 -1 1 0 0; 0 0 -4 0 0];
b_true = [3;2;5];
Ae_true = [1 0 0 0 0; 0 0 0 1 1];
be_true = [0;1];
P_true = polytope(A_true,b_true,Ae_true,be_true);
% check for equality
res(end+1,1) = isequal(P_,P_true,1e-10);


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
