function res = test_reachSet_shiftTime
% test_reachSet_shiftTime - unit test function for shiftTime
%
% Syntax:
%    res = test_reachSet_shiftTime()
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
% Written:       05-June-2023
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

% shift reach set by time
shift = 1;
R_shifted = shiftTime(R,shift);

% compare results
res = true;
for i=1:steps-1
    if ~withinTol(R_shifted.timePoint.time{i},...
            R.timePoint.time{i}+shift,1e-14)
        res = false; break
    end
end

% ------------------------------ END OF CODE ------------------------------
