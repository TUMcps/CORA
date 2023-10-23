function res = testMOSEK_polytope_compact
% testMOSEK_polytope_compact - unit test function of compact; comparing
%    computation using linprog from MOSEK and built-in MATLAB
%
% Syntax:
%    res = testMOSEK_polytope_compact
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
n = 5;
A = [eye(n);-eye(n);eye(n);-eye(n)];
b = [ones(2*n,1);2*ones(2*n,1)];
P{1} = polytope(A,b);

A = [2 1; -1 3; -2 -2; 1 -4];
b = ones(4,1);
P{2} = polytope(A,b);

% add redundant halfspace
A = [2 1; -1 3; -2 -2; 1 -4; 1 0];
b = ones(5,1);
P{3} = polytope(A,b);

% empty polytope: x1 + x2 >= 2, x1 - x2 >= 2, x1 <= -1
A = [-1 1; 0 1; 1 0; -1 -2; -1/sqrt(2) -1/sqrt(2)];
b = [0; 1; 2; 0; -sqrt(2)];
P{4} = polytope(A,b);

% check full-dimensionality using MOSEK linprog
P_mosek = cell(length(P),1);
for i=1:length(P)
    P_mosek{i,1} = compact(P{i});
end

% remove MOSEK from path
[suc,path2mosek] = removeSolverFromPath('mosek');
if ~suc
    res = false; return
end

% check full-dimensionality using MATLAB linprog
P_matlab = cell(length(P),1);
for i=1:length(P)
    P_matlab{i,1} = compact(P{i});
end

% add MOSEK to the path again
for i=1:length(path2mosek)
    addpath(path2mosek{i});
end

% compare results
res = true;
for i=1:length(P)
    if ~eq(P_mosek{i},P_matlab{i},1e-8)
        res = false; break
    end
end

% ------------------------------ END OF CODE ------------------------------
