function res = testMOSEK_polytope_isFullDim
% testMOSEK_polytope_isFullDim - unit test function of isFullDim; comparing
%    computation using linprog from MOSEK and built-in MATLAB
%
% Syntax:
%    res = testMOSEK_polytope_isFullDim
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

% Authors:       Mark Wetzlinger
% Written:       14-December-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% ensure that MOSEK is on the path
if ~isSolverInstalled('mosek')
    res = false; return
end

% init polytopes

% full-dimensional
A = [-1 -1; 1 0;-1 0; 0 1; 0 -1];
b = [2; 3; 2; 3; 2];
P{1} = polytope(A,b);

% degenerate: -2 <= x1 <= 2, x2 = 2
A = [1 0; -1 0; 0 1; 0 -1];
b = [2; 2; 2; -2];
P{2} = polytope(A,b);

% degenerate and unbounded
A = [1 1 0; 1 -1 0; -1 0 0];
b = zeros(3,1);
P{3} = polytope(A,b);

% check full-dimensionality using MOSEK linprog
res_mosek = false(length(P),1);
for i=1:length(P)
    res_mosek(i,1) = isFullDim(P{i});
end

% remove MOSEK from path
[suc,path2mosek] = removeSolverFromPath('mosek');
if ~suc
    res = false; return
end

% check full-dimensionality using MATLAB linprog
res_matlab = false(length(P),1);
for i=1:length(P)
    res_matlab(i,1) = isFullDim(P{i});
end

% add MOSEK to the path again
for i=1:length(path2mosek)
    addpath(path2mosek{i});
end

% compare results
res = all(res_mosek == res_matlab);

% ------------------------------ END OF CODE ------------------------------
