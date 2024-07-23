function res = test_polytope_projectHighDim
% test_polytope_projectHighDim - unit test function of projectHighDim
%
% Syntax:
%    res = test_polytope_projectHighDim
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

% Authors:       Tobias Ladner, Mark Wetzlinger
% Written:       18-September-2023
% Last update:   14-July-2024 (MW, vertex instantiation, special case)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true(0);

% 1D, fully empty
A = zeros(0,1); b = zeros(0,0);
P = polytope(A,b);
P_high = projectHighDim(P,3,2);
Ae_true = [1 0 0; 0 0 1]; be_true = [0;0];
P_true = polytope([],[],Ae_true,be_true);
res(end+1,1) = isequal(P_high,P_true);

% 2D, vertex instantiation
V = [2 1; -1 0; 0 -2]';
P = polytope(V);
P_high = projectHighDim(P,3,[3,1]);
V_true = [1 0 2; 0 0 -1; -2 0 0]';
V_high = vertices(P_high);
res(end+1,1) = compareMatrices(V_true,V_high);

% 2D, special case (no higher-dimensional space, just re-ordering)
A = [1 0; -1 1; -1 -1]; b = [1;1;1];
P = polytope(A,b);
P_high = projectHighDim(P,2,[2,1]);
A_true = [0 1; 1 -1; -1 -1]; b_true = [1;1;1];
P_true = polytope(A_true,b_true);
res(end+1,1) = isequal(P_high,P_true);

% 2D, both representations given
A = [1 0; -1 1; -1 -1]; b = [1;1;1];
P = polytope(A,b);
V = vertices(P);
P_high = projectHighDim(P,3,[1,3]);
V_high = vertices(P_high);
A_true = [1 0 0; -1 0 1; -1 0 -1]; b_true = [1;1;1];
Ae_true = [0 1 0]; be_true = 0;
V_true = [1 0 2; 1 0 -2; -1 0 0]';
P_true = polytope(A_true,b_true,Ae_true,be_true);
res(end+1,1) = isequal(P_high,P_true);
res(end+1,1) = compareMatrices(V_true,V_high);


% 2D, bounded
A = [1 0;-1 0;0 1;0 -1;1 1]; b = [1;1;1;1;1];
P = polytope(A,b);
P_high = projectHighDim(P,10,[4,5]);
% check support function for new dimensions
for i=[1 2 3 6 7 8 9 10]
    ei = zeros(10,1);
    ei(i) = 1;
    res(end+1) = supportFunc(P_high,ei) == 0;
    ei(i) = -1;
    res(end+1) = supportFunc(P_high,ei) == 0;
end
% check for projected dimensions
for i=[4,5]
    ei = zeros(10,1);
    ei(i) = 1;
    res(end+1) = supportFunc(P_high,ei) == 1;
    ei(i) = -1;
    res(end+1) = supportFunc(P_high,ei) == 1;
end

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
