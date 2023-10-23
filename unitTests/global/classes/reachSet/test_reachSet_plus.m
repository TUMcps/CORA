function res = test_reachSet_plus
% test_reachSet_plus - unit test function for plus
%
% Syntax:
%    res = test_reachSet_plus()
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

% instantiate set, propagation matrix
Z = zonotope([1;1],[1 0 -2; 2 -1 1]);
s = [1; -1];
dt = 0.02;
steps = 10;

% propagate sets
for i=1:steps
    timePoint.set{i,1} = Z + i*s;
    timePoint.time{i,1} = i*dt;
end
R = reachSet(timePoint);

% mapped reach set
R_shifted = R + s;

% compare results
res = true;
for i=1:steps-1
    if ~isequal(R.timePoint.set{i+1},R_shifted.timePoint.set{i},1e-14)
        res = false; break
    end
end

% ------------------------------ END OF CODE ------------------------------
