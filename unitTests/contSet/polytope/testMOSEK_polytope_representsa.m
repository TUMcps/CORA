function res = testMOSEK_polytope_representsa
% testMOSEK_polytope_representsa - unit test function of representsa;
%    comparing computation using linprog for emptiness check from MOSEK and
%    built-in MATLAB
%
% Syntax:
%    res = testMOSEK_polytope_representsa
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
 
% assume true
res = true;

% ensure that MOSEK is on the path
if ~isSolverInstalled('mosek')
    assert(false);
end

% init polytopes

% 1D: x <= 1, x >= 3
A = [1; -1];
b = [1; -3];
P{1} = polytope(A,b);

% polytope enclosing the origin
A = [2 1; -2 3; -2 -2; 4 1];
b = ones(4,1);
P{2} = polytope(A,b);

% empty polytope: x1 + x2 >= 2, x1 - x2 >= 2, x1 <= -1
A = [-1 -1; -1 1; 1 0];
b = [-2; -2; -1];
P{3} = polytope(A,b);

% check full-dimensionality using MOSEK linprog
res_mosek = false(length(P),1);
for i=1:length(P)
    res_mosek(i,1) = isFullDim(P{i});
end

% remove MOSEK from path
[suc,path2mosek] = removeSolverFromPath('mosek');
if ~suc
    assert(false);
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
assert(all(res_mosek == res_matlab));

% ------------------------------ END OF CODE ------------------------------
