function res = testLong_converter_commonocean2cora_NYM()
% testLong_converter_commonocean2cora_NYM - test for CommonOcean conversion
%
% Syntax:
%    res = testLong_converter_commonocean2cora_NYM()
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
    commonocean2cora('USA_NYM-1_20190613_T-1',false);


% expected results
lw = 1;
ls = 0;
ldo = 4776;
lso = 6;

x0_x = 581530.07;
x0_y = 4506018.3;
x0_time = 12;
x0_velocity = 0.15432;
x0_orientation = 1.859787;

goalSet_set_set_Vertices = [581532,4506004.2;581532,4506008.2;581547,4506008.2;581547,4506004.2];
goalSet_set_set_NumHoles = 0;
goalSet_set_set_NumRegions = 1;
goalSet_time = interval(12,346);
goalSet_orientation = interval(1.7859156,3.356712);

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
