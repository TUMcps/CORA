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
A = 1; b = 0;
P1 = polytope(A,b);
A = [1 1; -1 -1]; b = [1;1];
P2 = polytope(A,b);
P_cartProd = cartProd(P1,P2);
A_true = [1, 0, 0; 0, 1, 1; 0,-1,-1]; b_true = [0;1;1];
P_true = polytope(A_true,b_true);
% check for equality
res(end+1,1) = isequal(P_cartProd,P_true,1e-10);


% 1D and 2D (fullspace and bounded)
A = zeros(0,1); b = zeros(0,0);
P1 = polytope(A,b);
A = [1 0; -1 1; -1 -1]; b = [1;1;1];
P2 = polytope(A,b);
P_cartProd = cartProd(P1,P2);
A_true = [0 1 0; 0 -1 1; 0 -1 -1]; b_true = [1;1;1];
P_true = polytope(A_true,b_true);
res(end+1,1) = isequal(P_cartProd,P_true,1e-10);


% 2D and 2D (both only inequalities)
A = [-1 -1; 1 0;-1 0; 0 1; 0 -1]; b = [2; 3; 2; 3; 2];
P1 = polytope(A,b);
A = [1;-1]; b = [3; 2];
P2 = polytope(A,b);
P_cartProd = cartProd(P1,P2);
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
res(end+1,1) = isequal(P_cartProd,P_true,1e-10);


% 1D and 3D (both only equalities)
Ae = 5; be = 1;
P1 = polytope([],[],Ae,be);
Ae = [1 0 1; -1 0 0; 0 1 1]; be = [2;4;-3];
P2 = polytope([],[],Ae,be);
P_cartProd = cartProd(P1,P2);
Ae_true = [1 0 0 0; 0 1 0 1; 0 -1 0 0; 0 0 1 1];
be_true = [0.2;2;4;-3];
P_true = polytope([],[],Ae_true,be_true);
% check for equality
res(end+1,1) = isequal(P_cartProd,P_true,1e-10);


% 2D (only equalities) and 3D (only inequalities)
Ae = [1 0; 0 1]; be = [1;1];
P1 = polytope([],[],Ae,be);
A = [1 0 0; 0 1 0; 0 0 1; -1 -1 -1]; b = ones(4,1);
P2 = polytope(A,b);
P_cartProd = cartProd(P1,P2);
A_true = [0 0 1 0 0; 0 0 0 1 0; 0 0 0 0 1; 0 0 -1 -1 -1];
b_true = ones(4,1);
Ae_true = [1 0 0 0 0; 0 1 0 0 0];
be_true = [1;1];
P_true = polytope(A_true,b_true,Ae_true,be_true);
% check for equality
res(end+1,1) = isequal(P_cartProd,P_true,1e-10);


% 1D and 4D (both with equalities and inequalities)
A = [1;2]; b = [4;20]; Ae = 1; be = 0;
P1 = polytope(A,b,Ae,be); % with redundancies
A = [2 1 0 0; -1 1 0 0; 0 -4 0 0]; b = [3;2;5]; Ae = [0 0 1 1]; be = 1;
P2 = polytope(A,b,Ae,be);
P_cartProd = cartProd(P1,P2);
A_true = [0 2 1 0 0; 0 -1 1 0 0; 0 0 -4 0 0];
b_true = [3;2;5];
Ae_true = [1 0 0 0 0; 0 0 0 1 1];
be_true = [0;1];
P_true = polytope(A_true,b_true,Ae_true,be_true);
% check for equality
res(end+1,1) = isequal(P_cartProd,P_true,1e-10);


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
