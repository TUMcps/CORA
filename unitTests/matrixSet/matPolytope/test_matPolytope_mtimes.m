function res = test_matPolytope_mtimes
% test_matPolytope_mtimes - unit test function for mtimes
% 
% Syntax:
%    res = test_matPolytope_mtimes
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

% Authors:       Tobias Ladner
% Written:       02-May-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% scalar
V = []; V(:,:,1) = 1; V(:,:,2) = -2;
matP = matPolytope(V);
matPmtimes = matP * 2;
assert(compareMatrices(V * 2,matPmtimes.V));
matPmtimes = 2 * matP;
assert(compareMatrices(V * 2,matPmtimes.V));

% nx1 vector
V = []; V(:,:,1) = [0; 1]; V(:,:,2) = [1; -1]; V(:,:,3) = [-2; 0];
matPvec = matPolytope(V);
M = [1 2; 4 5; 5 6];
matPmtimes = M * matPvec;
assert(compareMatrices(pagemtimes(M,V),matPmtimes.V));

% matrix
V = [];
V(:,:,1) = [0 2; 1 -1; 1 -2];
V(:,:,2) = [1 1; -1 0; -2 1];
V(:,:,3) = [-2 0; 0 1; 1 -1];
matP = matPolytope(V);
M = [1 2 3; 4 5 6; 5 6 7];
matPmtimes = M * matP;
assert(compareMatrices(pagemtimes(M,V),matPmtimes.V));

% matP * matPvec
matPmatPvec = matP * matPvec;
assert(all(dim(matPmatPvec) == [3,1]));

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
