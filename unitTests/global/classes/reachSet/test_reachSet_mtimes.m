function res = test_reachSet_mtimes
% test_reachSet_mtimes - unit test function for mtimes
%
% Syntax:
%    res = test_reachSet_mtimes()
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
% Written:       10-November-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% init result
res = true;

% instantiate set, propagation matrix
Z = zonotope([1;1],[1 0 -2; 2 -1 1]);
A = [-1 -4; 4 -1];
dt = 0.02;
steps = 10;

% propagate sets
for i=1:steps
    timePoint.set{i,1} = expm(A*i*dt) * Z;
    timePoint.time{i,1} = i*dt;
end
R = reachSet(timePoint);

% mapped reach set
R_mapped = expm(A*dt) * R;

% compare results
for i=1:steps-1
    if ~isequal(R.timePoint.set{i+1},R_mapped.timePoint.set{i},1e-14)
        res = false; break
    end
end

% ------------------------------ END OF CODE ------------------------------
