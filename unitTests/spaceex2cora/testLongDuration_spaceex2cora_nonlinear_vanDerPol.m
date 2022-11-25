function pass = testLongDuration_spaceex2cora_nonlinear_vanDerPol()
% testLongDuration_spaceex2cora_nonlinear_vanDerPol - example of nonlinear reachability
%
% Syntax:  
%    testLongDuration_spaceex2cora_nonlinear_vanDerPol
%
% Inputs:
%    no
%
% Outputs:
%    pass - boolean 
% 
% Author:       Matthias Althoff, Raja Judeh
% Written:      26-June-2009
% Last update:  13-September-2018
% Last revision:15-January-2019 (SX_HA to SX_nonLin)

%------------- BEGIN CODE --------------


%% Specify continuous dynamics

vanDerPolSys = nonlinearSys(@vanderPolEq); %original system

%% Convert the system to xml file to test cora2spaceex and then back to CORA

cora2spaceex(vanDerPolSys,'vanDerPolSys');
spaceex2cora('vanDerPolSys.xml');
vanDerPolSys_SX = vanDerPolSys(); 

%% Simulation
R0 = zonotope([1.4 0.3 0; 2.3 0 0.05]);

params.x0 = randPoint(R0);
params.tFinal = 6.74;

% Simulating the original system
[simResOriginal.t, simResOriginal.x] = simulate(vanDerPolSys, params);
                             
% Simulating the converted system
[simResSX.t, simResSX.x] = simulate(vanDerPolSys_SX, params);

%% Compute error between the simulation of the two files

diff = simResOriginal.x - simResSX.x;
error = sqrt(sum(diff.^2,2));

if any(error > 1e-5) 
    disp('Failed Conversion: error = ' + string(max(error)));
    pass = false;
    return
end 

disp('Successful Conversion: error = ' + string(max(error)))
pass = true;

end

%------------- END OF CODE --------------