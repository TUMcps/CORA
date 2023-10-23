function res = testMOSEK_polytope_supportFunc
% testMOSEK_polytope_supportFunc - unit test function of supportFunc;
%    comparing computation using linprog from MOSEK and built-in MATLAB
%
% Syntax:
%    res = testMOSEK_polytope_supportFunc
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

% infeasible constraints
A = [-1 -1; -1 1; 1 0];
b = [-2; -2; -1];
P{1} = polytope(A,b);
dir{1} = [1;0];

% non-degenerate
A = [2 1; 1 3; 2 -2; -2 -2; -1 2];
A = (A' ./ vecnorm(A'))';
b = [2; 1; 2; 2; 2];
P{2} = polytope(A,b);
dir{2} = A';

% unbounded
A = [1 0; -1 0];
b = ones(2,1);
P{3} = polytope(A,b);
dir{3} = [0 1; 0 -1]';

% check full-dimensionality using MOSEK linprog
res_mosek = cell(length(P),1);
for i=1:length(P)
    for l=1:size(dir{i},2)
        res_mosek{i,1}(l) = supportFunc(P{i},dir{i}(:,l));
    end
end

% remove MOSEK from path
[suc,path2mosek] = removeSolverFromPath('mosek');
if ~suc
    res = false; return
end

% check full-dimensionality using MATLAB linprog
res_matlab = cell(length(P),1);
for i=1:length(P)
    for l=1:size(dir{i},2)
        res_matlab{i,1}(l) = supportFunc(P{i},dir{i}(:,l));
    end
end

% add MOSEK to the path again
for i=1:length(path2mosek)
    addpath(path2mosek{i});
end

% compare results
res = true;
for i=1:length(res_mosek)
    for l=1:length(res_mosek{i})
        if ( isfinite(res_mosek{i}(l)) && ~withinTol(res_mosek{i}(l),res_matlab{i}(l)) ) ...
                || ( ~isfinite(res_mosek{i}(l)) && res_mosek{i}(l) ~= res_matlab{i}(l) )s
            res = false; return
        end
    end
end

% ------------------------------ END OF CODE ------------------------------
