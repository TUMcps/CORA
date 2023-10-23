function res = test_polytope_mtimes
% test_polytope_mtimes - unit test function of linear map
%
% Syntax:
%    res = test_polytope_mtimes
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
% Written:       30-November-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = [];

% init polytope
V = [3 2; 0 3; -3 0; -1 -2; 2 -2]';
P = polytope(V);

% compute linear map
A = [2 1; -1 0];
P_mapped = A*P;

% compute vertices
V_mapped = vertices(P_mapped);

% map vertices
V_ = A*V;
P_ = polytope(V_);

% 1. check
res(end+1,1) = compareMatrices(V_,V_mapped,1e-14);

% 2. check: normalize and compare
temp1 = [P_mapped.A, P_mapped.b]';
temp1 = temp1 ./ vecnorm(temp1);
temp2 = [P_.A, P_.b]';
temp2 = temp2 ./ vecnorm(temp2);
res(end+1,1) = compareMatrices(temp1,temp2,1e-14);

% unbounded polytope
P = polytope([1 0; -1 0; 0 1], [1;1;1]);
P_ = A * P;
res(end+1,1) = ~isBounded(P_);

% degenerate polytope
P = polytope([1 0; -1 0; 0 1; 0 -1], [1;1;1;-1]);
P_ = A * P;
res(end+1,1) = ~isFullDim(P_);

% standard scaling matrix
P = polytope([1 0;-1 0; 0 1; 0 -1],[2;2;2;2]);
A = [2 0; 0 4];
P_ = mtimes(A, P);
P_true = polytope([1 0;-1 0; 0 1; 0 -1],[4;4;8;8]);
res(end+1,1) = P_ == P_true;


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
