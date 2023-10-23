function res = testLong_spaceex2cora_hybrid_bouncingball()
% testLong_spaceex2cora_hybrid_bouncingball - example for hybrid dynamics
%
% Syntax:
%    res = testLong_spaceex2cora_hybrid_bouncingball
%
% Inputs:
%    -
%
% Outputs:
%    res - boolean 

% Authors:       Matthias Althoff, Raja Judeh
% Written:       27-July-2016
% Last update:   13-September-2018
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%% Specify hybrid automation

% continuous dynamics 
A = [0 1; 0 0];
B = [0; 0];
c = [0; -9.81];
linSys = linearSys('linearSys',A,B,c);

% system parameters
alpha = -0.75;                  % rebound factor

% invariant set 
inv = polytope([-1,0],0);

% guard sets
guard = conHyperplane([1,0],0,[0,1],0);

% reset function
reset.A = [0, 0; 0, alpha]; reset.c = zeros(2,1);

% transitions
trans = transition(guard,reset,1);

% location object
loc = location('loc1',inv,trans,linSys); 

% hybrid automata
HA = hybridAutomaton(loc);

%% Convert the system to xml file to test cora2spaceex and then back to CORA

cora2spaceex(HA,'hybrid_bball');
spaceex2cora('hybrid_bball.xml');
HA_SX = hybrid_bball();


%% Simulate hybrid automations
R0 = zonotope([1;0],diag([0.05,0.05]));

params.x0 = randPoint(R0);
params.tFinal = 1.7;
params.startLoc = 1;

[~,x] = simulate(HA, params); %simulate the original file
[~,x_SX] = simulate(HA_SX, params); %simulate the converted file


%% Compute error between the simulation of the two files

% Note: the number of the channels in simResOriginal exceeds the number of
% channels in simResConverted by one (the first cell is the different one).
% Therefore, we avoid the first cell and compare the rest of the cells. 

for cell = 1:length(x)
    diff = x{cell} - x_SX{cell};
    error = vecnorm(diff); %has two values: one for each channel
    
    if any(error > 1e-5) 
        disp('Failed Conversion: error = ' + string(error));
        res = false;
        return
    end 
end

disp('Successful Conversion: error = ' + string(sum(error)));
res = true;

% ------------------------------ END OF CODE ------------------------------
