function res = testLong_spaceex2cora_hybrid_spacecraft()
% testLong_spaceex2cora_hybrid_spacecraft - test SpaceEx conversion for a
%                                       hybrid system with nonlinear guards
%
% Syntax:
%    res = testLong_spaceex2cora_hybrid_spacecraft()
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false 

% Authors:       Niklas Kochdumper
% Written:       18-May-2020
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Hybrid Automaton --------------------------------------------------------

HA = spacecraft_levelSet();


% SpaceEx Conversion ------------------------------------------------------

% generate SpaceEx model for hybrid automaton
cora2spaceex(HA,'spacecraft');

% convert the SpaceEx model to a CORA model
spaceex2cora('spacecraft.xml')
HA_SX = spacecraft(); 


% Simulation --------------------------------------------------------------

R0 = zonotope([[-900; -400; 0; 0],diag([25;25;0;0])]);

params.x0 = randPoint(R0);
params.tFinal = 200;
params.startLoc = 1;

% Simulating the original system
[~,x] = simulate(HA, params);
                             
% Simulating the converted system
[~,x_SX] = simulate(HA_SX, params);


% Compute Error -----------------------------------------------------------

for i = 1:length(x)
    
    diff = x{i} - x_SX{i};
    error = sum(diff.^2,1);
    
    if any(error > 1e-5) 
        disp('Failed Conversion: error = ' + string(max(error)));
        res = false;
        return
    end 
end

disp('Successful Conversion: error = ' + string(max(error)))
res = true;

% ------------------------------ END OF CODE ------------------------------
