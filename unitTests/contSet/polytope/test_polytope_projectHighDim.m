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

% Authors:       Tobias Ladner
% Written:       18-September-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true(0);

% 1D, fully empty
A = zeros(0,1); b = zeros(0,0);
P = polytope(A,b);
P_high = projectHighDim(P,3,2);
A_true = [1 0 0; 0 0 1]; b_true = [0;0];
P_true = polytope(A_true,b_true);
res(end+1,1) = isequal(P_high,P_true);


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
for i=[4 5]
    ei = zeros(10,1);
    ei(i) = 1;
    res(end+1) = supportFunc(P_high,ei) == 1;
    ei(i) = -1;
    res(end+1) = supportFunc(P_high,ei) == 1;
end

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
