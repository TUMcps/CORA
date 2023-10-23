function res = testLong_spaceex2cora_linear_building()
% testLong_spaceex2cora_linear_building - example of linear reachability 
% analysis from the ARCH17 friendly competition (building example)
%
% Syntax:
%    res = testLong_spaceex2cora_linear_building
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false 

% Authors:       Matthias Althoff, Raja Judeh
% Written:       09-February-2017
% Last update:   02-September-2018
% Last revision: 15-January-2019 (changed SX_HA to SX_lin)

% ------------------------------ BEGIN CODE -------------------------------

%% Specify continuous dynamics

G = load('build'); %loads the struct 'G' from file 'build'

% Initialize system
buildingSys = linearSys('buildingSys', G.A, G.B); %original system

%% Convert the system to xml file to test cora2spaceex and then back to CORA
cora2spaceex(buildingSys,'linear_building');
spaceex2cora('linear_building.xml');
buildingSysSX = linear_building(); 

%% Compute reachable set using zonotopes

R0 = interval([0.0002*ones(10,1); -0.0001; zeros(37,1)], ... 
              [0.00025*ones(10,1); 0.0001; zeros(37,1)]);
          
params.x0 = randPoint(R0);
params.u = 0.9;
params.tStart = 0; %start time
params.tFinal = 2; %final time

% Simulate the original system
[~,x] = simulate(buildingSys, params);

% Simulate the system converted from spaceEX
[~,xSX] = simulate(buildingSysSX, params);

%% Compute error between the simulation of the two files

num_channels = size(x, 2); %number of channels in the system

for ch = 1:num_channels
    diff = x(:,ch) - xSX(:,ch);
    error = norm(diff); %length of the difference vector

    if(error > 1e-5)
        disp('Failed Conversion: error = ' + string(error));
        res = false;
        return
    end
end

disp('Successful Conversion: error = ' + string(error))
res = true;

% ------------------------------ END OF CODE ------------------------------
