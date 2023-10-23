function res = testMOSEK_polytope_center
% testMOSEK_polytope_center - unit test function of center; comparing
%    computation using linprog from MOSEK and built-in MATLAB
%
% Syntax:
%    res = testMOSEK_polytope_center
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

% instantiate polytope
A = [1 1; -1 1; 1 -1; -1 -1];
b = [1; 1; 1; 1];
P = polytope(A,b);

% compute center using MOSEK linprog
c_mosek = center(P);

% remove MOSEK from path
[suc,path2mosek] = removeSolverFromPath('mosek');
if ~suc
    res = false; return
end

% compute center using MATLAB linprog
c_matlab = center(P);

% add MOSEK to the path again
for i=1:length(path2mosek)
    addpath(path2mosek{i});
end

% compare results
res = all(withinTol(c_matlab,c_mosek));

% ------------------------------ END OF CODE ------------------------------
