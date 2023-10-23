function res = test_reachSet_isequal
% test_reachSet_isequal - unit test function for isequal
%
% Syntax:
%    res = test_reachSet_isequal()
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
% Written:       01-May-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% instantiate set, propagation matrix
Z = zonotope([1;1],[1 0 -2; 2 -1 1]);
shift = [1e-6;-1e-8];
A = [-1 -4; 4 -1];
dt = 0.02;
dt_shift = 1e-10;
steps = 10;

% propagate sets
timePoint.set = cell(steps,1);
timePoint.time = cell(steps,1);
timePoint_.set = cell(steps,1);
timePoint_.time = cell(steps,1);
timePoint__.set = cell(steps,1);
timePoint__.time = cell(steps,1);

% time-point solution
for i=1:steps
    timePoint.set{i,1} = expm(A*i*dt) * Z;
    timePoint.time{i,1} = i*dt;
    % small shift (both time)
    timePoint_.set{i,1} = expm(A*i*dt) * Z + shift;
    timePoint_.time{i,1} = i*dt + dt_shift;
    % small shift (only time)
    timePoint__.set{i,1} = expm(A*i*dt) * Z;
    timePoint__.time{i,1} = i*dt + dt_shift;
end

timeInt.set = cell(steps-1,1);
timeInt.time = cell(steps-1,1);
timeInt_.set = cell(steps-1,1);
timeInt_.time = cell(steps-1,1);
timeInt__.set = cell(steps-1,1);
timeInt__.time = cell(steps-1,1);
% time-interval solution
for i=1:steps-1
    timeInt.set{i,1} = enclose(timePoint.set{i},timePoint.set{i+1});
    timeInt.time{i,1} = interval(timePoint.time{i},timePoint.time{i+1});
    % small shift (both)
    timeInt_.set{i,1} = enclose(timePoint_.set{i},timePoint_.set{i+1});
    timeInt_.time{i,1} = interval(timePoint_.time{i},timePoint_.time{i+1});
    % small shift (only time)
    timeInt__.set{i,1} = enclose(timePoint.set{i},timePoint.set{i+1});
    timeInt__.time{i,1} = interval(timePoint_.time{i},timePoint_.time{i+1});
end

% init reach set objects

% standard
R1 = reachSet(timePoint,timeInt);
% time-interval solution: set and time shifted
R2 = reachSet(timePoint,timeInt_);
% time-point solution: set and time shifted
R3 = reachSet(timePoint_,timeInt);
% both solutions: set and time shifted
R4 = reachSet(timePoint_,timeInt_);
% only time-point solution given
R5 = reachSet(timePoint);
% time-point solution: only time shifted
R6 = reachSet(timePoint__,timeInt);
% time-interval solution: only time shifted
R7 = reachSet(timePoint,timeInt__);

% compare reach set objects
res = isequal(R1,R1);
res(end+1,1) = ~isequal(R1,R2);
res(end+1,1) = ~isequal(R1,R2,1e-10);
res(end+1,1) = isequal(R1,R2,1e-4);
res(end+1,1) = ~isequal(R1,R3);
res(end+1,1) = isequal(R1,R3,1e-4);
res(end+1,1) = ~isequal(R1,R4);
res(end+1,1) = isequal(R1,R4,1e-4);
res(end+1,1) = ~isequal(R1,R5);
res(end+1,1) = ~isequal(R1,R6);
res(end+1,1) = isequal(R1,R6,1e-4);
res(end+1,1) = ~isequal(R1,R7);
res(end+1,1) = isequal(R1,R7,1e-4);

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
