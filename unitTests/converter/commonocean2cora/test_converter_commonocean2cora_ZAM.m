function res = test_converter_commonocean2cora_ZAM()
% test_converter_commonocean2cora_ZAM - test for CommonOcean conversion
%
% Syntax:
%    res = test_converter_commonocean2cora_ZAM()
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
    commonocean2cora('ZAM_Tutorial-2_1_T-1',false);


% expected results
lw = 1;
ls = 1;
ldo = 273;
lso = 3;

x0_x = 0;
x0_y = 0;
x0_time = 0;
x0_velocity = 15;
x0_orientation = 1.57079632;

goalSet_set_set_Vertices = [7.142857142857141,150;-7.142857142857146,150;-7.142857142857141,250;7.142857142857146,250];
goalSet_set_set_NumHoles = 0;
goalSet_set_set_NumRegions = 1;
goalSet_time = interval(-47,73);
goalSet_orientation = interval(1.27079632,1.87079632);


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
