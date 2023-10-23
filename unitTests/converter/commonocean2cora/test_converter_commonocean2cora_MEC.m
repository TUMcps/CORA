function res = test_converter_commonocean2cora_MEC()
% test_converter_commonocean2cora_MEC - test for CommonOcean conversion
%
% Syntax:
%    res = test_converter_commonocean2cora_MEC()
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false 

% Authors:       Bruno Maione
% Written:       14-March-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% CommonOcean conversion
[statObs,dynObs,x0,goalSet,waterways,shallows] = ...
    commonocean2cora('USA_MEC-1_20190129_T-22',false);


% expected results
lw = 0;
ls = 0;
ldo = 582;
lso = 0;

x0_x = 10234.782;
x0_y = -8690.8556;
x0_time = 21;
x0_velocity = 5.19544;
x0_orientation = 2.280253;

goalSet_set_set_Vertices = [2497.9801,681.14751;2497.9801,710.14751;2705.9801,710.14751;2705.9801,681.14751];
goalSet_set_set_NumHoles = 0;
goalSet_set_set_NumRegions = 1;
goalSet_time = interval(2610,2810);
goalSet_orientation = interval(2.13579589,2.3357959);


% check conversion
res = length(waterways) == lw;
res(end+1,1) = length(shallows) == ls;
res(end+1,1) = length(dynObs) == ldo;
res(end+1,1) = length(statObs) == lso;

res(end+1,1) = x0.x == x0_x && x0.y == x0_y && x0.time == x0_time ...
    && x0.velocity == x0_velocity && x0.orientation == x0_orientation;

res(end+1,1) = isequal(goalSet{1}.set.set.Vertices,goalSet_set_set_Vertices) ...
    && goalSet{1}.set.set.NumHoles == goalSet_set_set_NumHoles ...
    && goalSet{1}.set.set.NumRegions == goalSet_set_set_NumRegions ...
    && isequal(goalSet{1}.time,goalSet_time) ...
    && isequal(goalSet{1}.orientation,goalSet_orientation);

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
